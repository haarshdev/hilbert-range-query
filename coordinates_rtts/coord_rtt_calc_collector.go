package main

import (
	"bytes"
	"encoding/json"
	"log"
	"net/http"
	"os/exec"
	"strings"
	"sync"
	"time"
	"math"

	"github.com/hashicorp/serf/client"
	"github.com/hashicorp/serf/coordinate"
)

// NodeInfo represents information for each node
type NodeInfo struct {
	Name   string                 `json:"name"`
	Addr   string                 `json:"addr"`
	Port   int                    `json:"port"`
	Status string                 `json:"status"`
	Tags   map[string]string      `json:"tags"`
	Coord  *coordinate.Coordinate `json:"coordinate,omitempty"`
	RTTs   map[string]float64     `json:"rtts,omitempty"`
}

// ClusterStatus holds the complete cluster snapshot
type ClusterStatus struct {
	Timestamp string     `json:"timestamp"`
	Nodes     []NodeInfo `json:"nodes"`
}

var (
	cachedStatus []byte
	cacheLock    sync.RWMutex
)

// ------------------- Coordinate Operations -------------------

func collectAllCoordinates(serfClient *client.RPCClient) ([]NodeInfo, string) {
	members, err := serfClient.Members()
	if err != nil {
		log.Printf("ERROR: failed to get Serf members: %v", err)
		return nil, ""
	}
	log.Printf("Found %d Serf members", len(members))

	// Identify local node
	cmd := exec.Command("./serf", "info")
	out, err := cmd.Output()
	if err != nil {
		log.Printf("ERROR: failed to run 'serf info': %v", err)
		return nil, ""
	}

	localName := ""
	lines := strings.Split(string(out), "\n")
	for _, line := range lines {
		line = strings.TrimSpace(line)
		if strings.HasPrefix(line, "name =") {
			parts := strings.SplitN(line, "=", 2)
			if len(parts) == 2 {
				localName = strings.TrimSpace(parts[1])
				break
			}
		}
	}

	if localName == "" {
		log.Printf("ERROR: could not determine local node name from serf info output:\n%s", string(out))
		return nil, ""
	}
	log.Printf("Local node name: %s", localName)

	// Collect node info for all members
	nodes := make([]NodeInfo, 0, len(members))
	for _, member := range members {
		coord, err := serfClient.GetCoordinate(member.Name)
		if err != nil {
			log.Printf("WARNING: failed to get coordinate for node %s: %v", member.Name, err)
		}
		nodes = append(nodes, NodeInfo{
			Name:   member.Name,
			Addr:   member.Addr.String(),
			Port:   int(member.Port),
			Status: member.Status,
			Tags:   member.Tags,
			Coord:  coord,
		})
	}

	log.Printf("Collected %d nodes", len(nodes))
	return nodes, localName
}

// ------------------- RTT via Coordinates -------------------

func computeRTT(a, b *coordinate.Coordinate) float64 {
	if a == nil || b == nil {
		return -1
	}

	if len(a.Vec) != len(b.Vec) {
		return -1
	}

	sumsq := 0.0
	for i := 0; i < len(a.Vec); i++ {
		diff := a.Vec[i] - b.Vec[i]
		sumsq += diff * diff
	}

	rtt := math.Sqrt(sumsq) + a.Height + b.Height
	
	// Apply adjustment only if result stays positive
    adjusted := rtt + a.Adjustment + b.Adjustment
    if adjusted > 0.0 {
        rtt = adjusted
    }

	// RTT is returned in seconds by model — convert to milliseconds
	return rtt * 1000.0
}

func collectLocalRTTs(nodes []NodeInfo, localName string) {
	for i := range nodes {
		if nodes[i].Name != localName || nodes[i].Status != "alive" {
			continue
		}

		rtts := make(map[string]float64)
		for j := range nodes {
			if nodes[j].Name == localName || nodes[j].Status != "alive" {
				continue
			}
			rtts[nodes[j].Name] = computeRTT(nodes[i].Coord, nodes[j].Coord)
		}
		nodes[i].RTTs = rtts
	}
}

// ------------------- Cluster Status -------------------

func computeClusterStatus(serfClient *client.RPCClient) ClusterStatus {
	nodes, localName := collectAllCoordinates(serfClient)
	if nodes == nil || localName == "" {
		return ClusterStatus{}
	}

	collectLocalRTTs(nodes, localName)

	return ClusterStatus{
		Timestamp: time.Now().Format(time.RFC3339),
		Nodes:     nodes,
	}
}

// ------------------- Background Updater -------------------

func startBackgroundRTTUpdate(serfClient *client.RPCClient, interval time.Duration) {
	go func() {
		for {
			status := computeClusterStatus(serfClient)

			var buf bytes.Buffer
			if err := json.NewEncoder(&buf).Encode(status); err == nil {
				cacheLock.Lock()
				cachedStatus = buf.Bytes()
				cacheLock.Unlock()
			}

			time.Sleep(interval)
		}
	}()
}

// ------------------- HTTP Handler -------------------

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

// ------------------- Main -------------------

func main() {
	serfClient, err := client.ClientFromConfig(&client.Config{
		Addr: "127.0.0.1:7373",
	})
	if err != nil {
		log.Fatalf("Failed to connect to Serf RPC: %v", err)
	}
	defer serfClient.Close()

	startBackgroundRTTUpdate(serfClient, 30*time.Second)

	http.HandleFunc("/cluster-status", statusHandler())
	log.Fatal(http.ListenAndServe(":8080", nil))
}
