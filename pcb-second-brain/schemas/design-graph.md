# DesignIR — the in-memory design graph

Produced by the parsers, consumed by every checker, generator and report.
Pydantic models in `pcbbrain-core`. Serialises to JSON; large boards are
stored as Parquet tables for corpus statistics.

```
Design
  meta: {name, source_path, eda_tool, eda_version, imported_from, fetched_at, licence, owner: bool}
  sheets: [Sheet {path, title, parent}]
  components: [Component {
      refdes, value, part_id | null, mpn_raw, manufacturer_raw,
      footprint_ref, symbol_ref, sheet, dnp: bool,
      placement: {x_mm, y_mm, rot_deg, layer} | null,
      properties: {...}}]
  nets: [Net {
      name, class: power|gnd|signal|diffpair|clock|analog|unknown,
      diffpair_partner | null, netclass_name, width_hint_mm | null}]
  connections: [{refdes, pin, net}]
  stackup: {layers: [{name, type, material, thickness_mm, cu_oz}], total_mm,
            impedance_table: [{geometry, width_mm, gap_mm, z0_ohm}]} | null
  board: {outline_polygons, size_mm, mounting_holes, keepouts} | null
  tracks: [{net, layer, width_mm, length_mm, start, end}]
  vias: [{net, x, y, drill_mm, size_mm, layers}]
  zones: [{net, layer, area_mm2, is_plane}]
  rules_in_use: {dru_path, netclasses: [...]}
  provenance: {parser, parser_version, parsed_at, sha256}
```

Derived views (computed lazily, cached per design):
- `pin_net_graph()` → NetworkX graph of pins and nets (bipartite).
- `neighbourhood(refdes, hops)` → subgraph used for motif mining and T3.
- `distance(refdes_a, refdes_b)` → placement distance, if placed.
- `net_stats(net)` → total length, min width, via count, layers, return-plane coverage.
- `power_tree()` → regulators and the rails they feed, inferred from part roles and nets.
