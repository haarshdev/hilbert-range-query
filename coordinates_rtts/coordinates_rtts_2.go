package main

import (
	"bytes"
	"encoding/json"
	"log"
	"math"
	"net/http"
	"sync"
	"time"
	"fmt"

	"github.com/hashicorp/serf/client"
	"github.com/hashicorp/serf/coordinate"
)

type Float3 float64

func (f Float3) MarshalJSON() ([]byte, error) {
    v := float64(f)
    if math.IsNaN(v) || math.IsInf(v, 0) {
        return []byte("null"), nil
    }
    // JSON NUMBER token with exactly 3 decimals (no quotes)
    return []byte(fmt.Sprintf("%.3f", v)), nil
}

type NodeInfo struct {
	Name   string                 `json:"name"`
	Addr   string                 `json:"addr"`
	Port   int                    `json:"port"`
	Status string                 `json:"status"`
	Tags   map[string]string      `json:"tags"`
	Coord  *coordinate.Coordinate `json:"coordinate,omitempty"`
	RTTs   map[string]Float3      `json:"rtts,omitempty"`
}

type ClusterStatus struct {
	Timestamp string     `json:"timestamp"`
	Nodes     []NodeInfo `json:"nodes"`
}

var (
	cachedStatus []byte
	cacheLock    sync.RWMutex
)

// dist computes the RTT between two coordinates using the Serf coordinate system.
// Returns RTT in milliseconds (float64).
func dist(a, b *coordinate.Coordinate) float64 {
	if a == nil || b == nil {
		return -1
	}

	if len(a.Vec) != len(b.Vec) {
		log.Printf("⚠️ Coordinate dimensionality mismatch between %v and %v", a, b)
		return -1
	}

	sumsq := 0.0
	for i := 0; i < len(a.Vec); i++ {
		diff := a.Vec[i] - b.Vec[i]
		sumsq += diff * diff
	}
	rtt := math.Sqrt(sumsq) + a.Height + b.Height

	// Apply adjustment and avoid negatives
	adjusted := rtt + a.Adjustment + b.Adjustment
	if adjusted > 0.0 {
		rtt = adjusted
	}

	// RTT is returned in seconds by model — convert to milliseconds
	return rtt * 1000.0
}

func truncate3(x float64) float64 {
	return math.Floor(x*1000) / 1000
}

func computeRTTData(serfClient *client.RPCClient) ClusterStatus {
	members, err := serfClient.Members()
	if err != nil {
		log.Printf("❌ Error retrieving members: %v", err)
		return ClusterStatus{}
	}

	nodes := make([]NodeInfo, 0, len(members))
	for _, member := range members {
		coord, _ := serfClient.GetCoordinate(member.Name)
		nodes = append(nodes, NodeInfo{
			Name:   member.Name,
			Addr:   member.Addr.String(),
			Port:   int(member.Port),
			Status: member.Status,
			Tags:   member.Tags,
			Coord:  coord,
		})
	}

	for i := range nodes {
		if nodes[i].Status != "alive" || nodes[i].Coord == nil {
			continue
		}

		rtts := make(map[string]Float3)
		for j := range nodes {
			if i == j || nodes[j].Status != "alive" || nodes[j].Coord == nil {
				continue
			}

			rttMs := dist(nodes[i].Coord, nodes[j].Coord)
			rtts[nodes[j].Name] = Float3(truncate3(rttMs))
		}
		nodes[i].RTTs = rtts
	}

	return ClusterStatus{
		Timestamp: time.Now().Format(time.RFC3339),
		Nodes:     nodes,
	}
}

func startBackgroundRTTUpdate(serfClient *client.RPCClient, interval time.Duration) {
	go func() {
		for {
			log.Println("🔄 RTT monitor active... computing coordinate-based RTTs")
			status := computeRTTData(serfClient)

			var buf bytes.Buffer
			err := json.NewEncoder(&buf).Encode(status)
			if err != nil {
				log.Printf("❌ Failed to encode status JSON: %v", err)
			} else {
				cacheLock.Lock()
				cachedStatus = buf.Bytes()
				cacheLock.Unlock()
				log.Printf("✅ Cluster status cached at %s with %d nodes", status.Timestamp, len(status.Nodes))
			}

			time.Sleep(interval)
		}
	}()
}

func statusHandler() http.HandlerFunc {
	return func(w http.ResponseWriter, r *http.Request) {
		cacheLock.RLock()
		defer cacheLock.RUnlock()

		if cachedStatus == nil {
			http.Error(w, "Status not ready yet", http.StatusServiceUnavailable)
			return
		}

		w.Header().Set("Content-Type", "application/json")
		w.Write(cachedStatus)
	}
}

func main() {
	serfClient, err := client.ClientFromConfig(&client.Config{
		Addr: "127.0.0.1:7373",
	})
	if err != nil {
		log.Fatalf("❌ Failed to connect to Serf RPC: %v", err)
	}
	defer serfClient.Close()

	startBackgroundRTTUpdate(serfClient, 30*time.Second)

	http.HandleFunc("/cluster-status", statusHandler())
	log.Println("🚀 Server running on :8080 (coordinate-based RTT mode)")
	log.Fatal(http.ListenAndServe(":8080", nil))
}
