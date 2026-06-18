# Cerebro FS — object durability & replication

Cerebro stores its read sets in a SeaweedFS cluster. By default the tiered
single-server models (`single-server`, `distributed-hpc`) run **one** volume
server with replication `000` — a single copy. A disk failure there loses data.
cerebro-fs provides a durable single-host posture and documents the multi-host upgrade.

## Single-host replication (`single-server-replicated`)

```bash
cerebro stack deploy \
  --fs-model single-server-replicated \
  --fs-hot /mnt/ssd0/hot   --fs-cold /mnt/hdd0/cold \
  --fs-replica-hot /mnt/ssd1/hot --fs-replica-cold /mnt/hdd1/cold \
  ...
```

This renders two tiered volume servers — `cerebro-fs-primary` and
`cerebro-fs-replica` — in the **same** data centre and rack, and sets replication
to `001` (one replica on a different server in the same rack). SeaweedFS then keeps
a second copy of every volume across the two disk sets, so the loss of one
server's disks does not lose data. The existing durability gate
(`validate_durability`) confirms `001` is satisfiable by the two servers before the
stack is written.

**Use physically separate disks.** `--fs-replica-hot`/`--fs-replica-cold` must be
on different physical devices from `--fs-hot`/`--fs-cold`. If both copies live on
the same disk they share a failure domain and the redundancy is only nominal.

**What this does and does not protect against.** It protects against **disk**
failure (one device dies, the other still holds a full copy). It does **not**
protect against **host** failure — both volume servers, the master, and the filer
run on one machine. For host-failure durability, see *Multi-host* below.

### SeaweedFS replication codes (recap)

A code is `xyz`: `x` copies on other data centres, `y` on other racks (same DC),
`z` on other servers (same rack). Total copies = `x + y + z + 1`.

| code | copies | placement | use |
|------|--------|-----------|-----|
| `000` | 1 | single server | no redundancy (default tiered) |
| `001` | 2 | +1 server, same rack | **single-host replica** |
| `010` | 2 | +1 rack | two racks, one site |
| `100` | 2 | +1 data centre | two sites / hosts labelled as DCs |

The durability gate refuses any code that needs more copies than the topology has
volume servers (SeaweedFS would otherwise reject writes).

## Monitoring under-replication

Replication is only as good as its observability — a silently under-replicated
volume is a latent single copy. Watch for it:

- **Master metrics** (Prometheus, `cerebro-fs-master:9324`) and per-volume-server
  metrics (`:9325` primary, `:9328` replica) expose volume counts and disk usage.
  Alert when the replica server's volume count diverges from the primary's, or when
  either server is unreachable.
- **`weed shell`** against the master gives the authoritative view:
  - `volume.list` — shows each volume's replica placement and whether it is
    under-replicated.
  - `volume.fix.replication` — re-replicates under-replicated volumes onto a
    server that can satisfy the code (used once a failed disk/server is replaced).
- A volume that shows fewer replicas than the code requires, or a volume server
  that has been down past your alert window, is the signal to investigate and, once
  hardware is restored, run `volume.fix.replication`. SeaweedFS re-replicates
  automatically when a target server can satisfy the code; the DR runbook (DR-1)
  documents the manual nudge and the reconcile-based verification that the copy
  healed.

## Multi-host (host-failure durability)

Single-host replication shares one machine. To survive losing a host, place the
replica volume server (and its disks) on a **second machine** and label the
topology so SeaweedFS spreads copies across the fault boundary:

1. Run `cerebro-fs-replica` on host B, reachable by the master on host A over the
   cluster network (open the volume server's port and the master gRPC port `19333`
   between hosts; keep them off the public internet).
2. Give the two servers labels that reflect the real boundary and raise the code to
   match:
   - different **rack**, same site → `010`;
   - different **data centre** / site → `100`.
   Set the matching `-dataCenter` / `-rack` on each volume server and
   `-defaultReplication` on the master.
3. The master and filer remain single points of failure until they too are made
   redundant — that is an HA concern beyond single-host replication (a SeaweedFS master quorum + a
   MongoDB replica set for the filer/catalogue), noted here as the next step and
   covered for the catalogue by the catalogue backups and the disaster-recovery runbook.

The deploy CLI currently renders the single-host layout; the multi-host layout is a
small generalisation of the same `ReplicaConfig` (distinct `center`/`rack` per
server) plus per-host compose files, planned as a follow-on once the single-host
path is validated on real hardware.
