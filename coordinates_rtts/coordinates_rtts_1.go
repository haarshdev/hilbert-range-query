package main

import (
	"bytes"
	"encoding/json"
	"fmt"
	"log"
	"net/http"
	"os/exec"
	"strconv"
	"strings"
	"sync"
	"time"

	"github.com/hashicorp/serf/client"
	"github.com/hashicorp/serf/coordinate"
)

type NodeInfo struct {
	Name   string                 `json:"name"`
	Addr   string                 `json:"addr"`
	Port   int                    `json:"port"`
	Status string                 `json:"status"`
	Tags   map[string]string      `json:"tags"`
	Coord  *coordinate.Coordinate `json:"coordinate,omitempty"`
	RTTs   map[string]float64     `json:"rtts,omitempty"`
}

type ClusterStatus struct {
	Timestamp string     `json:"timestamp"`
	Nodes     []NodeInfo `json:"nodes"`
}

var (
	cachedStatus []byte
	cacheLock    sync.RWMutex
)

func getRealRTT(from, to string) (float64, error) {
	cmd := exec.Command("./serf", "rtt", from, to)
	output, err := cmd.CombinedOutput()

	if err != nil {
		log.Printf("RTT command failed: %v, output: %s", err, output)
		return -1, fmt.Errorf("serf rtt command failed: %s", string(output))
	}

	parts := strings.Fields(string(output))
	if len(parts) < 6 {
		return -1, fmt.Errorf("unexpected output: %s", output)
	}

	rttMs, err := strconv.ParseFloat(parts[5], 64)
	if err != nil {
		return -1, fmt.Errorf("invalid RTT value: %s", parts[5])
	}

	return rttMs, nil
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
		if nodes[i].Status != "alive" {
			continue
		}

		rtts := make(map[string]float64)
		for j := range nodes {
			if i == j || nodes[j].Status != "alive" {
				continue
			}

			rttMs, err := getRealRTT(nodes[i].Name, nodes[j].Name)
			if err != nil {
				rtts[nodes[j].Name] = -1
				continue
			}

			rtts[nodes[j].Name] = rttMs
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
			log.Println("🔄 RTT monitor active... computing RTTs")
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
	log.Println("🚀 Server running on :8080")
	log.Fatal(http.ListenAndServe(":8080", nil))
}
