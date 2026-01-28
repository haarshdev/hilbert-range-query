package main

import (
	"bytes"
	"encoding/json"
	"log"
	"math"
	"net"
	"net/http"
	"sync"
	"time"

	"github.com/hashicorp/serf/client"
	"github.com/hashicorp/serf/coordinate"
)

// ------------------- Data Model -------------------

type RTTStatus struct {
	Timestamp string             `json:"timestamp"`
	Node      string             `json:"node"`
	Addr      string             `json:"addr"`
	RTTs      map[string]float64 `json:"rtts"`
}

var (
	cachedStatus []byte
	cacheLock    sync.RWMutex
)

// ------------------- RTT Computation -------------------

func computeRTT(a, b *coordinate.Coordinate) float64 {
	if a == nil || b == nil || len(a.Vec) != len(b.Vec) {
		return -1
	}
	sumsq := 0.0
	for i := range a.Vec {
		diff := a.Vec[i] - b.Vec[i]
		sumsq += diff * diff
	}
	rtt := math.Sqrt(sumsq) + a.Height + b.Height
	adj := rtt + a.Adjustment + b.Adjustment
	if adj > 0 {
		rtt = adj
	}
	return rtt * 1000.0 // seconds → ms
}

// ------------------- Coordinate Fetch -------------------

func collectLocalRTTs(serfClient *client.RPCClient) RTTStatus {
	members, err := serfClient.Members()
	if err != nil {
		log.Printf("ERROR: get members: %v", err)
		return RTTStatus{}
	}

	stats, err := serfClient.Stats()
	if err != nil {
		log.Printf("ERROR: get stats: %v", err)
		return RTTStatus{}
	}
	localName := stats["agent"]["name"]

	localMember := findMemberByName(members, localName)
	localCoord, err := serfClient.GetCoordinate(localName)
	if err != nil {
		log.Printf("ERROR: get local coordinate: %v", err)
		return RTTStatus{}
	}

	rtts := make(map[string]float64)
	for _, m := range members {
		if m.Status != "alive" || m.Name == localName {
			continue
		}
		coord, err := serfClient.GetCoordinate(m.Name)
		if err != nil {
			continue
		}
		rtt := computeRTT(localCoord, coord)
		// Round to 3 decimals
		rtts[m.Name] = math.Round(rtt*1000) / 1000
	}

	return RTTStatus{
		Timestamp: time.Now().Format(time.RFC3339),
		Node:      localName,
		Addr:      localMember.Addr.String(),
		RTTs:      rtts,
	}
}

func findMemberByName(members []client.Member, name string) client.Member {
	for _, m := range members {
		if m.Name == name {
			return m
		}
	}
	// Fallback empty member if not found
	return client.Member{Addr: net.ParseIP("0.0.0.0")}
}

// ------------------- Background Updater -------------------

func startBackgroundRTTUpdate(serfClient *client.RPCClient, interval time.Duration) {
	go func() {
		for {
			status := collectLocalRTTs(serfClient)

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
			http.Error(w, "RTT data not ready", http.StatusServiceUnavailable)
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

	startBackgroundRTTUpdate(serfClient, 10*time.Second)

	http.HandleFunc("/rtts", statusHandler())
	log.Println("Serving RTT data on :8080/rtts")
	log.Fatal(http.ListenAndServe(":8080", nil))
}
