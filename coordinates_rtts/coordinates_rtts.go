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

var (
	cachedStatus []byte
	cacheLock    sync.RWMutex
)

// getRealRTT executes the CLI command `./serf rtt A B` and extracts the RTT in ms
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

// computeRTTData builds the full NodeInfo slice with RTT values
func computeRTTData(serfClient *client.RPCClient) []NodeInfo {
	members, err := serfClient.Members()
	if err != nil {
		log.Printf("❌ Error retrieving members: %v", err)
		return nil
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

	return nodes
}

// startBackgroundRTTUpdate runs periodically to recompute full cluster RTTs
func startBackgroundRTTUpdate(serfClient *client.RPCClient, interval time.Duration) {
	go func() {
		for {
			log.Println("🔄 Refreshing cluster RTT data...")
			nodes := computeRTTData(serfClient)

			var buf bytes.Buffer
			err := json.NewEncoder(&buf).Encode(nodes)
			if err != nil {
				log.Printf("❌ Failed to encode status JSON: %v", err)
			} else {
				cacheLock.Lock()
				cachedStatus = buf.Bytes()
				cacheLock.Unlock()
				log.Printf("✅ Cluster status cached with %d nodes", len(nodes))
			}

			time.Sleep(interval)
		}
	}()
}

// statusHandler returns the latest cached JSON instantly
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
