# Rule taxonomy (W2)

The kernel's canonical design-rule taxonomy: one set of rule kinds that is the union of KiCad's custom rule
language (`.kicad_dru`) and board-setup constraints, Altium Designer's design-rule types and Eagle's design rules,
with a mapping in both directions. Rule records use `schemas/rule.schema.json`; the objects they talk about are
those of `schemas/design-graph.md` (DesignIR).

Status: draft for K1 (docs/09-roadmap.md, "taxonomy covers 100% of rule types in golden boards"). Written
2026-11-04 in work item W2, task S1; verifier corrections applied the same day in task F1.

## 1. Purpose, sources and how to read this file

### 1.1 Purpose

- Give every rule the kernel checks a stable canonical id (section 3), so imported KiCad, Altium and Eagle rules,
  fab capability floors (W7) and interface rules (W8) all land on the same kinds.
- Define one scope language for "which objects does this rule apply to" (section 2.4) and how KiCad conditions,
  Altium queries and Eagle classes translate into it (section 2.5).
- List every native rule item found in the W2 source notes with its canonical target (section 4), so an importer
  or exporter can be written from the tables without re-reading tool sources.
- Show that every rule and constraint kind used by the ten golden boards maps to a canonical kind (section 6).

### 1.2 Sources

| # | Source | What it gave | Limits |
|---|---|---|---|
| S1 | KiCad source, GitLab commit `4c23271608ae7fa9e584c6f63b57e54163832c48` (master, 2026-11-03), sparse checkout of 258 rule-related files (191 under `pcbnew/`, chiefly `pcbnew/drc/` and the Altium and Eagle importers; 30 under `common/`; 21 `qa/`; 11 `include/`; 4 `eeschema/erc/`; 1 `tools/`) | DRU keywords, parser, engine resolution order, implicit rules, board-setup keys, netclasses, DRC providers and codes, ERC codes and pin matrix, KiCad's Altium and Eagle importers | Files outside the sparse paths were fetched read-only from GitLab at the same commit where needed; they are marked with a dagger (†). The legacy KiCad 5 board and `.pro` parsers are not in the checkout. |
| S2 | Rule-bearing settings extracted from the ten golden boards (`golden-rules/<board>.txt`) and, where the extract truncated or filtered, the boards' own files at their pinned commits (`golden/<board>/src/...`) | What real boards use: KiCad `.kicad_pro` rules and netclasses, KiCad 5 `(setup)` and `(net_class)`, KiCad 5 `.pro`, Eagle `<designrules>` and `.dru`, the 37 Altium `Rules6` records of odrive-v3 and its `.PrjPcb` ERC settings | The extractor (scratchpad script `scripts/golden_rules.py`, not committed) misses unnamed Eagle `<designrules>` blocks (its line 54) and cuts long records; the board files were read for those cases. |
| S3 | Public third-party files on GitHub and crates.io (keys W1-W13, table 1.4) | Real Altium rule records (a `.RUL` export, `Rules6` dumps, ASCII PcbDocs), Altium's own DelphiScript sample with the `TRuleKind` names, third-party format notes | unverified (L1: confirm via Altium scripting, Q-A12) and W13 are someone else's reverse engineering; see Q-A12 in section 8. |
| S4 | WebSearch result summaries of Altium documentation pages (keys AD1-AD15) | Category lists, UI names, units of rules with no record in any file, Room Definition semantics | altium.com is blocked from this environment; no page was opened. Every AD-keyed fact needs checking in L1. |
| S5 | One third-party Eagle glossary (key OMX) | Meanings of Eagle `<param>` names | No official Eagle help could be read (the Debian help mirror, forum.kicad.info and forums.autodesk.com failed DNS). |

Blocked or missing: docs.kicad.org and altium.com (network policy); official Eagle documentation; a running KiCad,
Altium or Eagle (none in the cloud phase). Nothing here was checked by running a tool; statements derived from
reading code are marked "inferred".

The working notes this file was built from (A1 KiCad DRU language, A2 KiCad setup and ERC, A3 Altium, A4 Eagle,
A5 golden usage) and the generator that checks every native item has a mapping live in the session scratchpad and
are not committed. Every citation below points at the primary source, so this file stands alone.

### 1.3 How to read the tables

- **Canonical id**: snake_case, stable. Kinds are grouped in fourteen categories (section 3): the thirteen of the W2
  brief, plus "General" for custom assertions, data-integrity checks and engine settings that fit none of them.
- **Target syntax** in the mapping tables: `kind` or `kind.param`, `kind[variant]` for kinds that group several
  native checks, and a bracket for origin and scope, for example `hole_size.min [board_floor, type(via.micro)]`.
  `@field` means a rule-model field (section 2.1) rather than a kind, for example `@severity`.
- **Quality**: `exact` = round trip without loss; `lossy: <why>` = expressible with the stated loss; `none` = the tool
  has no equivalent (a "none" may name an approximation the exporter can use). `model` = the item maps onto a
  rule-model field, not a kind.
- **Arity**: `U` unary (one object), `B` binary (a pair of objects), `-` not a checked rule (generation or preference).
- **Units** are explicit everywhere: mm, mil, deg, %, count, ratio (fraction, 0.25 = 25 %), ohm, V, s, ps, mm2.
  Conversions: 1 mil = 0.0254 mm; KiCad internal unit 1 nm; Altium internal unit 0.0001 mil = 2.54 nm (common/io/altium/altium_props_utils.cpp:113, mils × 10000).
  The DRU parser multiplies only `hole_aspect_ratio` by 1000 (pcbnew/drc/drc_rule_parser.cpp:612, :922-923); a DRU
  `solder_paste_rel_margin` value goes through the length unit resolver and is not scaled (:614-619, :921-924). The
  paste ratio is held multiplied by 1000 only on the engine's item-override and board-default paths
  (pcbnew/drc/drc_engine.cpp:1396, :2081), which is the mismatch behind Q-K1.

### 1.4 Provenance forms and keys

| Form | Meaning |
|---|---|
| `pcbnew/...:N`, `common/...:N`, `eeschema/...:N`, `include/...:N`, `qa/...:N`, `tools/...:N` | Line N of the KiCad checkout (S1). A trailing † marks a file fetched from GitLab at the same commit. |
| `erc_item.cpp#LN` | websearch: https://gitlab.com/kicad/code/kicad/-/blob/4c23271608ae7fa9e584c6f63b57e54163832c48/eeschema/erc/erc_item.cpp#LN (fetched at the same commit; not in the checkout) |
| `golden-rules/<board>.txt:N` | Line N of the golden-board extract (S2). |
| `golden/<board>/src/<path>[:N]` (also tagged "golden-src") | The board's own file at its pinned commit, relative to `pcb-second-brain/`. Used only where the extract cut, filtered or skipped data. This is a fourth form beyond the three in the brief; it is a plain file:line. |
| `golden/odrive-v3/src/v3/PCB.PcbDoc rec <NAME>` | The `Rules6/Data` record named NAME, re-read without the extract's 400-character cut. |
| `odrive-DRC Dn` | Item Dn of the `Design Rule Checker Options6/Data` stream of the same PcbDoc (decoded list in section 4.2.4). |
| `Wn:L`, `Wn rec X` | websearch: URL of Wn (table below), line L, or record X in a CR-terminated ASCII PcbDoc. |
| `ADn` | websearch: URL of ADn (Altium documentation; known from search summaries only). |
| `OMX` | websearch: https://github.com/Omniblox/Eagle-Loader/blob/master/eagle/eagle-xml-notes.md (third-party Eagle glossary written against EAGLE 7.5.0 Light; fetched 2026-11-04). |
| `[E]`, `[X] raw :N`, `[P] raw :N`, `[M]`, `[M3]`, `[M5]`, `[D]` | Eagle value sources: [E] the two EEZ "Mechanical" boards (golden-rules/eez-psu-h24005.txt:48, :52); [X] `golden/eez-psu-h24005/src/Consolidated/Eagle files/EEZ PSU consolidated r5B13a.brd` line N; [P] the EEZ panel `.brd` line N; [M] moteus controller boards (golden-rules/moteus.txt:91); [M3] r4.3 (:103); [M5] r4.5 (:107); [D] `hw/controller/MF_4_mil_drill_4_Layer.dru` (golden-rules/moteus.txt:115). |
| `unverified (GK)` | unverified (general knowledge; verify in L1 against the local Altium install) |
| `unverified (GK-E)` | unverified (general knowledge; verify in L1 against a local EAGLE install) |
| `GK-K` | unverified (general knowledge; the KiCad legacy parser is not in the checkout) |
| `pcb-second-brain/docs/...` | Project requirement documents (used only for fab-only kinds that no tool defines). |
| `Q-xn` | An open question in section 8. |

Web keys (all URLs pinned to a commit where the host allows it):

| Key | URL | What it is |
|---|---|---|
| W1 | https://github.com/Altium-Designer-addons/scripting-reference/blob/a29001c80ae2d5fe6a17213e3f0ee7b6e77616c8/Delphiscript%20Scripts/PCB/PCB%20Rules/IterateRules.pas | Altium Limited DelphiScript sample (2008): `TRuleKind` 0-51 to `RULEKIND` strings, lines 16-67 |
| W5 | https://github.com/florentcournoyer16/S7_PCB/blob/cca53fcd4ec97899f231150f984f4589c8e32b90/01%20-%20SCHEMAS%20ET%20PCB/99%20-%20WORKING%20COPY/S7CAS-PLEIADES-CARTE-MERE_P01/S7CAS-PLEIADES-CARTE-MERE_P20_Design-Rules.RUL | Altium `.RUL` rule export, 44 records |
| W6 | https://github.com/erik-edmonds/openbms/blob/af6e8cde9f6e5e65be79de02f06790f775c2e7b6/hardware/worker/Sources/slave-PcbDoc/Root%20Entry/Rules6/Data.dat.txt | `Rules6/Data` text dump |
| W7 | https://github.com/nihospr01/OpenSpeechPlatform-2020B/blob/592fd4da3741f976be8293419e986093473fb95f/Hardware/DJB%20(Digital%20BTE-RIC)/djb_v2/djb_v2_flex/Outputs/Assembly/djb_v2_flex_ASCII.PcbDoc | ASCII PcbDoc (query and legacy scope forms) |
| W8 | https://github.com/novena-next/docs/blob/d8d495bd5c94d431d2a6c79dd2584008e7ae8e05/design_source/lidsense_board/ASCII-laptop-lidsense1.PcbDoc | ASCII PcbDoc with a Room Definition |
| W9 | https://github.com/alejandromb/octopus-drone/blob/113838e3a7db6e9346c2c34888fe4a0f41b1b20f/hardware/ESP32_S2_Drone_V1_2/src/PCB_Mainboard_ESP32_S2_Drone_V1_2.PcbDoc | ASCII PcbDoc with schematic-generated Supply Nets rules |
| W10 | https://github.com/BogdanSorlea/GSOC/blob/63b1aecf1024213b4efb81ded85742a6b377b388/Layout%20Files/PCB_Project/test/board-PcbDoc/Root%20Entry/Rules6/Data.dat.txt | `Rules6` dump with a Room Definition |
| W11 | https://github.com/iryek-rpi/pico-encoder/blob/64ff03c7bc49d9a10345be3cee1c845e9a8f2747/hw/parts/W5500-EVB-Pico/altium2kicad/W5500-EVB-Pico_V100-PcbDoc/Root%20Entry/Rules6/Data.dat.txt | `Rules6` dump (Board Outline Clearance) |
| W12 | https://github.com/HOPE-MICROELECTRONICS/Product_Portfolio/blob/804c2a081acae20833192178a1b7295ae91fadf8/Sub-1GHz/CMT2300A/SCH_PCB/CMT2310A.pcbdoc | ASCII PcbDoc (kind names only) |
| W13 | https://github.com/lgili/Fenolite/blob/b2b0b06347adf77c477b27e5da1768ddb2a2791c/docs/formats/altium/pcb-copper.md | Another project's Altium format notes (mostly labelled INFERRED by its author) |
| AD1 | https://www.altium.com/documentation/altium-designer/pcb/design-rule-types | "65 types of design rules, in 10 categories" |
| AD2 | https://altium.com/documentation/node/248066/printable/print?version=21 | Category lists incl. High Speed and Manufacturing members |
| AD3 | https://altium.com/documentation/altium-designer/pcb-design-rule-types?version=17 | SMT, Mask, Plane, Testpoint categories |
| AD4 | https://www.altium.com/documentation/15.1/display/ADES/PCB_Dlg-Form_DesignRules((PCB+Rules+and+Constraints+Editor))_AD | UI names of Electrical, Routing, High Speed, Placement members |
| AD5 | https://www.altium.com/cn/documentation/altium-designer/pcb-electrical-rule-creepage-distance?version=21 | Creepage Distance |
| AD6 | https://altium.com/cn/documentation/altium-designer/pcb-high-speed-rule-return-path | Return Path (unary) |
| AD7 | https://altium.com/documentation/altium-designer/pcb-high-speed-rule-max-via-stub-length-back-drilling | Max Via Stub Length (Back Drilling) |
| AD8 | https://altium.com/documentation/node/248641/printable/print ; https://www.altium.com/documentation/altium-designer/pcb-high-speed-rules | Parallel Segment, Matched Lengths, Daisy Chain Stub Length |
| AD9 | https://www.altium.com/documentation/altium-designer/pcb/design-rule-types/signal-integrity | Signal Integrity rules and units |
| AD10 | https://altium.com/documentation/altium-designer/pcb-smt-rules | SMD To Corner, SMD To Plane, SMD Neck-Down, SMD Entry |
| AD11 | https://www.altium.com/documentation/altium-designer/pcb-placement-rule-height ; https://www.altium.com/documentation/knowledge-base/altium-designer/how-to-create-an-area-on-the-pcb-where-component-footprint-placement-is-restricted-based-upon-footprint-height | Height, Room Definition |
| AD12 | https://www.altium.com/vi/altium-designer/whats-new ; https://www.altium.com/jp/altium-designer/coming-soon | Z-Axis Clearance (AD26), wire bonding |
| AD13 | https://techdocs.altium.com/display/ADOH/Scoping+Design+Rules ; https://www.eevblog.com/forum/altium/altium-query-to-find-free-pads/ | Scope is a query; priority 1 is highest; first matching rule wins |
| AD14 | https://altium.com/documentation/altium-designer/project-options-connection-matrix?version=19 ; https://altium.com/documentation/altium-designer/project-options-error-reporting?version=19 | ERC connection matrix and error levels |
| AD15 | https://www.altium.com/documentation/altium-designer/pcb-placement-rule-room-definition?version=21 ; https://altium.com/documentation/altium-designer/pcb-placement-rule-room-definition?version=17 | Room Definition: a component scoped by several Room Definition rules obeys all of them (search summary, 2026-11-04) |

## 2. Canonical rule model

### 2.1 What a rule is

A canonical rule is: **kind** + **scope** (unary, or binary with a pair relation) + **layer filter** +
**parameters** (min / opt / max with a unit, or a typed value) + **severity** + **priority** + **origin** +
**provenance**. The record format is `schemas/rule.schema.json`; the fields added for this model are optional, so
existing W8-style rule records stay valid (section 2.8).

| Field | Meaning | Values | Notes |
|---|---|---|---|
| `kind` | Canonical id from section 3 | snake_case | Absent on rules that bind only to a `check` callable (pure interface rules). |
| `variant` | Sub-check of a grouping kind | snake_case, listed per kind in section 3 | For example `hole_in_courtyard` + `npth`. |
| `scope.expr` | Canonical scope expression (section 2.4) for the first (or only) object | string | Binary rules: either `expr` over A and `expr_b` over B (Altium form), or one `expr` that names `A.` and `B.` (KiCad form). |
| `scope.expr_b` | Scope of the second object of a binary rule | string | Forbidden for unary kinds. |
| `scope.pair` | Relations the pair must satisfy | `same_net`, `different_net`, `coupled_diffpair`, `same_component`, `different_component`, `same_layer`, `adjacent_layer` | Altium `NETSCOPE` and `LAYERKIND` land here. |
| `scope.symmetric` | Match the pair in either order | bool, default true | KiCad evaluates (A, B) then (B, A) (pcbnew/drc/drc_rule_condition.cpp:105-117). |
| existing `scope.*` selectors | `part`, `part_family`, `package`, `interface`, `net_class`, `fab`, `process`, `std`, `stage` | as before | All given selectors and `expr` must match. `net_class: x` equals `expr: netclass('x')`. |
| `layers` | Layer filter | canonical layer ids or tokens (section 2.3) | Omitted = all layers. KiCad `(layer ...)`; Altium per-layer parameter families; Eagle layer groups. |
| `constraint` | Parameters | `min`, `opt`, `max` (numbers) with `unit` (`mm`, `deg`, `ps`, `s`, `count`, `ratio`, `percent`, `ohm`, `V`, `mm2`, `none`); `value` for tokens, flags, lists or expressions; `options` for named switches | `opt` is the target or preferred value (router default, tuning target). |
| `severity` | Report level | `blocker`, `error`, `warning`, `info`, plus `ignore` and `inherit` | Section 2.7. |
| `priority` | Order among rules of one kind and origin | integer >= 1; 1 is highest | Altium `PRIORITY` copied; KiCad file order reversed (section 5). |
| `origin` | Which tier the rule belongs to | `fab_floor`, `board_floor`, `tool_default`, `netclass`, `item_default`, `explicit`, `item_override` | Section 5. `item_default` is an object-level setting that applies only when no explicit rule matches (KiCad footprint zone connection, a zone's own connection, thermal gap and spoke width; section 5, P9). |
| `floor_exempt` | The rule may go below the floors | bool, default false | Set on imported KiCad DRU rules so results match KiCad (section 5, P8). |
| `enabled` | Rule is active | bool, default true | Altium `ENABLED`; Altium DRC check lists; Eagle `check*` switches. |
| `native` | Where the rule came from | `tool`, `kind`, `name`, `text`, `version` | Kept for round trips; never interpreted by checks. |
| `check`, `params`, `provenance`, `confidence`, `tags`, `title`, `description` | as before | | `check` is `pcbbrain_checks.drc.<kind>` or `pcbbrain_checks.erc.<kind>` by convention. Tool-specific extras with no canonical field go in `params`. Tag `waive` marks a rule whose findings are emitted as dismissed (KiCad severity `exclusion`). |

### 2.2 Arity and pair relations

- **Unary** kinds test one object (track width, hole size, text height, via count of a net). Their scope uses only A.
  KiCad rejects a condition naming B when any constraint of the rule is unary, except `AB.isCoupledDiffPair()` for
  length, net-chain length and skew (pcbnew/drc/drc_rule_parser.cpp:385-406; pcbnew/drc/drc_rule.cpp:88-118); the
  canonical model has the same rule. KiCad's unary list is `DRC_CONSTRAINT::IsUnary()` (pcbnew/drc/drc_rule.cpp:88-118);
  a constraint outside it accepts B even when it measures one geometric property.
- **Binary** kinds test a pair (clearances, diff-pair gap, courtyard overlap, zone connection of a pad to a zone).
  `connection_width` is binary with an optional B: KiCad evaluates it with the first two items that contribute copper
  at the narrow point, for example a pad and a zone, and B is absent when only one item contributes
  (pcbnew/drc/drc_test_provider_connection_width.cpp:450-454).
  Altium resolves rules with a one-object and a two-object call (unverified (L1: confirm via Altium scripting, Q-A12)); the two scope queries are
  `SCOPE1EXPRESSION` and `SCOPE2EXPRESSION`. Only Clearance (binary, W6:1) and Return Path (unary, AD6) have a sourced
  arity on the Altium side; every other Altium arity is unverified (GK) and follows the canonical kind.
- **Pair relations** are separate from the two scopes: `different_net` (Altium `DifferentNets`, KiCad `clearance`
  semantics), `same_net` (Altium `SameNetOnly`), `coupled_diffpair` (KiCad `AB.isCoupledDiffPair()`), `same_layer` and
  `adjacent_layer` (Altium `LAYERKIND`), `same_component` (Altium's pad-to-pad-in-footprint exception, KiCad net ties).

### 2.3 Objects, layers and units

Object types are those of DesignIR (`schemas/design-graph.md`); types marked (+) are not in that file yet and
must be added in W6 (docs/01-architecture.md already lists pins, footprints and pads in the model).

| Canonical type (subtypes) | DesignIR source | KiCad `A.Type` label / class | Altium | Eagle |
|---|---|---|---|---|
| `board` | `board` | (board) | `Board` scope | board |
| `net`, `netclass`, `diffpair` | `nets[]`, `Net.netclass_name`, `Net.diffpair_partner` | `NetName`, `NetClass`; diff pairs from `_P`/`_N` names | `InNet`, `InNetClass`, `InDifferentialPair` (named pair objects) | `signal`, `class` |
| `net_chain` (+) | (+) | net chains (master) | none | none |
| `component`, `footprint` | `components[]`, `Component.footprint_ref` | 'Footprint' / `FOOTPRINT` | `InComponent`, `HasFootprint` | element / package |
| `courtyard` (+) (top, bottom) | (+) | F.CrtYd / B.CrtYd shapes | component bodies (3D) | tKeepout / bKeepout (GK-E) |
| `pad` (+) (`smd`, `pth`, `npth`, `connector`, `aperture`; fabrication role `bga`, `fiducial`, `testpoint`, `heatsink`, `castellated`, `pressfit`) | `connections[]` (refdes, pin) | 'Pad' / `PAD`; `Pad_Type`: Through-hole, SMD, Edge connector, NPTH mechanical, Non-component; `Fabrication_Property` (pcbnew/pad.cpp:3649†, :3978†) | `IsPad`; pad classes | pad (THT), smd |
| `hole` (+) (`plated`, `unplated`, `slot`) | via `drill_mm` only | holes of pads and vias | Hole object | drill / hole |
| `track` (`segment`, `arc`) | `tracks[]` | 'Track' (segments, arcs and teardrop zones) | Track, Arc | wire |
| `via` (`through`, `blind`, `buried`, `micro`; `stitching`) | `vias[]` | 'Via'; `Via_Type`: Through, Blind, Buried, Micro (pcbnew/pcb_track.cpp:3684†) | `IsVia`, `IsStitchingVia` | via |
| `zone` (`copper`, `plane`, `fill`) | `zones[]` (`is_plane`) | 'Zone' | Polygon, plane layers, Region | polygon |
| `keepout` / `area` (+) (rule area, room, restrict shape) | `board.keepouts` | rule areas | Room, keepout primitives | tRestrict, bRestrict, vRestrict |
| `text` (+) (`field`, `box`), `graphic` (+), `dimension` (+), `barcode` (+), `group` (+), `table` (+), `generator` (+) | (+) | 'Text', 'Text Box', 'Graphic', 'Dimension', 'Leader', 'Barcode', 'Group', 'Table', 'Generator' (pcbnew/pcbexpr_evaluator.cpp:780-788) | Text, Fill | text, wire on non-signal layers |
| `edge` | `board.outline_polygons` | Edge.Cuts and Margin | board outline | Dimension layer (20) |
| `sheet`, `symbol`, `pin`, `label`, `wire`, `bus` | `sheets[]`, `connections[]`; label, wire, bus (+) | schematic objects (ERC) | schematic objects | schematic objects |

**Layers.** Canonical layer ids are positional, so they survive tool-specific names: `cu:1` (top copper) to
`cu:N` (bottom copper), with tokens `outer` (`cu:1`, `cu:N`), `inner` (`cu:2`..`cu:N-1`), `copper` (all) and `all`;
non-copper layers `silk:top|bottom`, `mask:top|bottom`, `paste:top|bottom`, `adhesive:top|bottom`,
`courtyard:top|bottom`, `fab:top|bottom`, `edge`, `margin`, `user:<name>`. Glob wildcards are allowed.
Translation: KiCad `F.Cu` -> `cu:1`, `In<k>.Cu` -> `cu:<k+1>`, `B.Cu` -> `cu:N`, `F.SilkS` -> `silk:top`,
`F.Mask` -> `mask:top`, `F.Paste` -> `paste:top`, `F.CrtYd` -> `courtyard:top`, `F.Fab` -> `fab:top`,
`Edge.Cuts` -> `edge`, `(layer outer)` / `(layer inner)` -> the tokens (pcbnew/drc/drc_rule_parser.cpp:929-978).
Altium `TOPLAYER` / `Top Layer` -> `cu:1`, `BOTTOMLAYER` -> `cu:N`, `MIDLAYER<k>` -> the stack position of that mid
layer, which needs the `Board6` layer stack (not yet decoded, Q-A11); `Top Overlay` -> `silk:top`; `OnMultiLayer` ->
through-hole objects. Eagle 1 -> `cu:1`, 16 -> `cu:N`, 2-15 -> inner layers in `layerSetup` order; restrict layers
41-43 (common/io/eagle/eagle_parser.h:1240-1242); other Eagle layer numbers unverified (GK-E).

**Units.** Lengths are stored in mm, angles in deg, delays in ps, SI times in s, counts unitless, ratios as
fractions, percentages as %, impedance in ohm, voltage in V, areas in mm2. Importers convert (section 1.3).

### 2.4 Canonical scope expression

A scope expression is a boolean expression over the item(s) under test. Grammar (EBNF):

```
expr        = or_expr ;
or_expr     = and_expr , { "or" , and_expr } ;
and_expr    = not_expr , { "and" , not_expr } ;
not_expr    = [ "not" ] , primary ;
primary     = "(" , expr , ")" | "all" | "none" | predicate | comparison ;
predicate   = [ item , "." ] , name , "(" , [ arg , { "," , arg } ] , ")" ;
arg         = string | quantity | item | name ;          (* name: a mode or side token, e.g. enclosed, front *)
comparison  = operand , cmp , operand ;
operand     = [ item , "." ] , value_fn , "(" , [ string ] , ")" | literal | operand , arith , operand ;
value_fn    = "prop" | "field" | "custom_property" | "pin_count" ;   (* value-returning functions, table below *)
item        = "A" | "B" | "A.parent" | "B.parent" ;
cmp         = "==" | "!=" | "<" | "<=" | ">" | ">=" | "~=" ;       (* "~=" is a glob match *)
arith       = "+" | "-" | "*" | "/" ;
literal     = quantity | string | "true" | "false" | "null" ;
quantity    = number , [ unit ] ;   unit = "mm" | "mil" | "in" | "um" | "deg" | "ps" | "fs" | "ohm" | "V" | "s" | "%" ;
string      = "'" , { ( char - "'" ) | "\'" } , "'" ;   (* "\'" is a quote inside the string; case-insensitive;
                                           "*" and "?" are wildcards in predicate arguments and on the right
                                           of "==" / "~=" *)
```

Precedence, from tightest: comparison, `not`, `and`, `or`. `not` applies to the whole `primary` that follows it, so
`not A.prop('Pad Type') == 'SMD'` means not(A.prop('Pad Type') == 'SMD'). KiCad differs in two places, so exporters
write explicit parentheses: (1) KiCad's `!` binds tighter than the comparison operators
(common/libeval_compiler/grammar.lemon:31-32: `%nonassoc` comparisons, then `%right G_BOOL_NOT`; in Lemon a later
declaration binds tighter), so
`!A.NetName == 'GND'` parses as `(!A.NetName) == 'GND'`; the KiCad exporter writes `!( ... )` around every negated term
that is not a single function call or identifier; (2) KiCad's grammar appears to bind `||` tighter than `&&`
(inferred from grammar.lemon:28-29), so the exporter parenthesises every mixed and/or (Q-K6).

Strings: the canonical form escapes a quote as `\'`, as KiCad does (common/libeval_compiler/libeval_compiler.cpp:224-229).
Export restrictions: KiCad cannot write a name that ends in a backslash (its escaper refuses it,
common/libeval_compiler/libeval_compiler.cpp:1472-1480; netclasses with such names get no implicit rules,
pcbnew/drc/drc_engine.cpp:270-274), and Altium queries have no escape at all (KiCad's tokenizer of Altium queries notes
"TODO: escaping?", pcbnew/pcb_io/altium/altium_rule_transformer.cpp:118), so a name containing `'` is not exported to
Altium; both cases are dropped and reported.

Predicates (default item A; `B.` selects the second object):

| Predicate | True when | KiCad | Altium | Eagle |
|---|---|---|---|---|
| `type(t, ...)` | object type or subtype is one of t (section 2.3), e.g. `type(via.micro)`, `type(pad.smd)` | `A.Type`, `A.Pad_Type`, `A.Via_Type`, `isMicroVia()` ... | `IsPad`, `IsVia`, object kinds | object kind of md* rows |
| `net(p)`, `netclass(p)`, `netclass_exact(p)`, `diffpair(p)` | net name / effective class contains / class equals / pair base name matches p | `NetName`, `hasNetclass`, `hasExactNetclass`, `inDiffPair` | `InNet`, `InNetClass`, `InDifferentialPair` | `signal/@class` |
| `net_chain(p)`, `net_chain_class(p)` | net belongs to a matching chain / chain class | `inNetChain`, `hasNetChain`, `inNetChainClass` | none | none |
| `component(p)`, `component_class(p)`, `footprint(p)`, `in_footprint(sel)`, `package_style(s)`, `pad_class(p)` | refdes / component class / library footprint / parent footprint / package style / pad class | `Reference`, `hasComponentClass`, `Library_Link`, `memberOfFootprint`; package style and pad classes: none | `InComponent`, `InComponentClass`, `HasFootprint`, `IsBGA` ..., `InPadClass` | element name |
| `group(p)`, `sheet(p[, children])`, `has_custom_property(k)` | grouping and metadata (boolean) | `memberOfGroup`, `memberOfSheet[OrChildren]`, `hasCustomProperty` | none | none |
| `on_layer(l)`, `layer_under_test(l)` | object exists on l / the layer being checked is l | `existsOnLayer`, `L` | `OnLayer`, `OnTopLayer` ... | layer number |
| `in_area(name, mode)` with mode `intersects`, `enclosed`, `intersects_fill` | object lies in a rule area, room or keepout | `intersectsArea`, `enclosedByArea`, `intersectsKeepout` | `WithinRoom`, `TouchesRoom` | restrict shapes |
| `in_courtyard(sel, side)` with side `any`, `front`, `back` | object touches the courtyard of a footprint | `intersects[Front\|Back]Courtyard` | none (bodies) | none |
| `from_to(a, b)` | object lies on the copper path between two pads | `fromTo` | `FromTo` scope kind | none |
| `plated()`, `stacked_via()` | plated hole / stacked micro-via | `isPlated`, `isStackedVia` | none | none |
| pair relations (`scope.pair`) | see section 2.2 | `AB.isCoupledDiffPair()`, `A.NetName == B.NetName` | `NETSCOPE`, `LAYERKIND` | class pairs |

Value functions (operands of a comparison, never booleans on their own):

| Function | Value | KiCad | Altium | Eagle |
|---|---|---|---|---|
| `prop('Display Name')` | any registered KiCad property, compared with units; the namespace is KiCad's property display names (Appendix A) | `A.<Name>` | object properties where a KiCad property has the same meaning (Altium pad `Name` -> `prop('Pad Number')`) | none |
| `field('name')` | text of a footprint field (string; empty when the object has no such field) | `A.getField('name')`, which pushes the field text (pcbnew/pcbexpr_functions.cpp:1629-1691; help: pcbnew/dialogs/panel_setup_rules_help_8expression_functions.md:40-43) | none | none |
| `custom_property('k')` | value of a custom property (string) | `A.customProperty('k')`, which pushes the property value (pcbnew/pcbexpr_functions.cpp:1823-1864; master only) | none | none |
| `pin_count()` | number of pads of the component (count) | none: KiCad registers no pad-count property; the exporter expands to an explicit refdes list | `CompPinCount` | none |

Semantics, chosen to match KiCad where it is explicit: string comparison is case-insensitive and a right-hand
string with `*` or `?` is a wildcard (common/libeval_compiler/libeval_compiler.cpp:132-137, :224-226); a property
missing on an object makes both `==` and `!=` false (:120-149), and the evaluator first retries the property on the
parent footprint (pcbnew/pcbexpr_evaluator.cpp:477-503); a binary rule matches the pair in either order unless
`scope.symmetric` is false; `${VAR}` text variables are expanded on import and the original text is kept in
`native.text` (pcbnew/drc/drc_engine.cpp:693-709).

### 2.5 Translating KiCad conditions, Altium queries and Eagle classes

**KiCad condition elements** (items, literals, operators, special fields and all 34 registered functions,
pcbnew/pcbexpr_functions.cpp:1921-1965):

| # | Kind | KiCad element | Canonical scope expression | Quality | Provenance |
|---|---|---|---|---|---|
| 1 | item | `A` | A | exact | pcbnew/dialogs/panel_setup_rules_help_3items.md:3 |
| 2 | item | `B` | B (binary rules: expr_b, or B. inside expr) | exact | pcbnew/dialogs/panel_setup_rules_help_3items.md:4-8 |
| 3 | item | `AB` | pair(...) predicates | exact | pcbnew/pcbexpr_evaluator.cpp:733-737 |
| 4 | item | `L` | @layers, or layer_under_test(l) | exact, except for edge_clearance and silk-to-edge, which KiCad evaluates with no layer under test, so a `(layer ...)` filter does not narrow them (lossy, section 5, P4) | pcbnew/dialogs/panel_setup_rules_help_3items.md:5; pcbnew/drc/drc_test_provider_edge_clearance.cpp:176, :218; pcbnew/drc/drc_engine.cpp:1819-1820 |
| 5 | literal | `null` | null | exact | pcbnew/pcbexpr_evaluator.cpp:698-703 |
| 6 | navigation | `A.Parent` / `B.Parent` / `A.Parent.<Field>` | A.parent / B.parent | exact | pcbnew/pcbexpr_evaluator.cpp:705-727 |
| 7 | literal | string `'...'` | '...' string (case-insensitive; * ? glob on the right operand) | exact | common/libeval_compiler/libeval_compiler.cpp:543-547 |
| 8 | literal | number with unit | quantity with unit (mm, mil, in, deg, ps, fs) | exact | pcbnew/pcbexpr_evaluator.cpp:857-899 |
| 9 | semantics | undefined values | three-valued: comparison with a missing property is false for == and != | exact | common/libeval_compiler/libeval_compiler.cpp:120-149 |
| 10 | operator | `\|\|` | or | lossy: KiCad \|\| binds tighter than && (inferred); importer must keep KiCad's parse tree | common/libeval_compiler/libeval_compiler.cpp:581-583 |
| 11 | operator | `&&` | and | lossy: see row 10 | common/libeval_compiler/libeval_compiler.cpp:576-578 |
| 12 | operator | `!` | not | lossy: KiCad's `!` binds tighter than comparisons (`!A.NetName == 'GND'` is `(!A.NetName) == 'GND'`); import keeps KiCad's parse tree, export writes `!( ... )` around every negated non-atomic term (section 2.4) | common/libeval_compiler/libeval_compiler.cpp:592; common/libeval_compiler/grammar.lemon:31-32, :139 |
| 13 | operator | `==` | == | exact | common/libeval_compiler/libeval_compiler.cpp:556-559 |
| 14 | operator | `!=` | != | exact | common/libeval_compiler/libeval_compiler.cpp:561-564 |
| 15 | operator | `<` `<=` `>` `>=` | < <= > >= | exact | common/libeval_compiler/libeval_compiler.cpp:566-574 |
| 16 | operator | `+` `-` `*` `/` | + - * / | exact | common/libeval_compiler/libeval_compiler.cpp:590-595 |
| 17 | operator | `( )` | ( ) | exact | common/libeval_compiler/libeval_compiler.cpp:598-599 |
| 18 | operator | `.` | . (member access) | exact | common/libeval_compiler/libeval_compiler.cpp:601 |
| 19 | operator | `,` | , (argument separator) | exact | common/libeval_compiler/libeval_compiler.cpp:602 |
| 20 | operator precedence | precedence (low to high) | canonical precedence comparison > not > and > or; exporter parenthesises every mixed and/or and every negated comparison | lossy: see rows 10-12 | common/libeval_compiler/grammar.lemon:26-36 |
| 21 | special field | `A.NetClass` / `B.NetClass` | netclass('x') (effective class contains x) | exact | pcbnew/pcbexpr_evaluator.cpp:753-761 |
| 22 | special field | `A.ComponentClass` / `B.ComponentClass` | component_class('x') | exact | pcbnew/pcbexpr_evaluator.cpp:762-770 |
| 23 | special field | `A.NetName` / `B.NetName` | net('x'); A.net == B.net -> pair same_net | exact | pcbnew/pcbexpr_evaluator.cpp:771-779 |
| 24 | special field | `A.Type` / `B.Type` | type(...) (label table in section 2.3) | exact | pcbnew/pcbexpr_evaluator.cpp:780-788 |
| 25 | semantics | property name resolution | prop('Display Name') | exact | pcbnew/pcbexpr_evaluator.cpp:804-841 |
| 26 | semantics | value conversion | value conversion rules of prop() | exact | pcbnew/pcbexpr_evaluator.cpp:532-616 |
| 27 | semantics | function names | function names case-insensitive | exact | pcbnew/pcbexpr_evaluator.cpp:679-689 |
| 28 | function | `A.existsOnLayer(...)` | on_layer('l') | exact | pcbnew/pcbexpr_functions.cpp:1921 |
| 29 | function | `A.isPlated()` | plated() | exact | pcbnew/pcbexpr_functions.cpp:1923 |
| 30 | function (deprecated) | `A.insideCourtyard(...)` | in_courtyard(sel, any) | exact | pcbnew/pcbexpr_functions.cpp:1927 |
| 31 | function (deprecated) | `A.insideFrontCourtyard(...)` | in_courtyard(sel, front) | exact | pcbnew/pcbexpr_functions.cpp:1928 |
| 32 | function (deprecated) | `A.insideBackCourtyard(...)` | in_courtyard(sel, back) | exact | pcbnew/pcbexpr_functions.cpp:1929 |
| 33 | function | `A.intersectsCourtyard(...)` | in_courtyard(sel, any) | exact | pcbnew/pcbexpr_functions.cpp:1930 |
| 34 | function | `A.intersectsFrontCourtyard(...)` | in_courtyard(sel, front) | exact | pcbnew/pcbexpr_functions.cpp:1931 |
| 35 | function | `A.intersectsBackCourtyard(...)` | in_courtyard(sel, back) | exact | pcbnew/pcbexpr_functions.cpp:1932 |
| 36 | function (deprecated) | `A.insideArea(...)` | in_area(name, intersects_fill) | exact | pcbnew/pcbexpr_functions.cpp:1934 |
| 37 | function | `A.intersectsArea(...)` | in_area(name, intersects) | exact | pcbnew/pcbexpr_functions.cpp:1935 |
| 38 | function | `A.intersectsKeepout(...)` | in_area(name, intersects_fill) | exact | pcbnew/pcbexpr_functions.cpp:1936 |
| 39 | function | `A.enclosedByArea(...)` | in_area(name, enclosed) | exact | pcbnew/pcbexpr_functions.cpp:1937 |
| 40 | function | `A.isMicroVia()` | type(via.micro) | exact | pcbnew/pcbexpr_functions.cpp:1939 |
| 41 | function | `A.isBlindVia()` | type(via.blind) | exact | pcbnew/pcbexpr_functions.cpp:1940 |
| 42 | function | `A.isBuriedVia()` | type(via.buried) | exact | pcbnew/pcbexpr_functions.cpp:1941 |
| 43 | function | `A.isBlindBuriedVia()` | type(via.blind) or type(via.buried) | exact | pcbnew/pcbexpr_functions.cpp:1942 |
| 44 | function | `A.isStackedVia()` | stacked_via() | exact | pcbnew/pcbexpr_functions.cpp:1943 |
| 45 | function (deprecated) | `A.memberOf(...)` | group('x') | exact | pcbnew/pcbexpr_functions.cpp:1945 |
| 46 | function | `A.memberOfGroup(...)` | group('x') | exact | pcbnew/pcbexpr_functions.cpp:1946 |
| 47 | function | `A.memberOfFootprint(...)` | in_footprint(sel) | exact | pcbnew/pcbexpr_functions.cpp:1947 |
| 48 | function | `A.memberOfSheet(...)` | sheet('p') | exact | pcbnew/pcbexpr_functions.cpp:1948 |
| 49 | function | `A.memberOfSheetOrChildren(...)` | sheet('p', children) | exact | pcbnew/pcbexpr_functions.cpp:1949 |
| 50 | function | `A.fromTo(...)` | from_to('REF-PAD', 'REF-PAD') | exact | pcbnew/pcbexpr_functions.cpp:1951 |
| 51 | function | `AB.isCoupledDiffPair()` | pair coupled_diffpair | exact | pcbnew/pcbexpr_functions.cpp:1952 |
| 52 | function | `A.inDiffPair(...)` | diffpair('base') | exact | pcbnew/pcbexpr_functions.cpp:1953 |
| 53 | function | `A.inNetChain(...)` | net_chain('x') | exact | pcbnew/pcbexpr_functions.cpp:1954 |
| 54 | function | `A.hasNetChain()` | net_chain('*') | exact | pcbnew/pcbexpr_functions.cpp:1955 |
| 55 | function | `A.inNetChainClass(...)` | net_chain_class('x') | exact | pcbnew/pcbexpr_functions.cpp:1956 |
| 56 | function | `A.getField(...)` | field('name') (value function, used inside a comparison: `A.getField('MPN') == 'X'` -> `field('MPN') == 'X'`) | exact | pcbnew/pcbexpr_functions.cpp:1958, :1629-1691 |
| 57 | function | `A.hasNetclass(...)` | netclass('x') | exact | pcbnew/pcbexpr_functions.cpp:1960 |
| 58 | function | `A.hasExactNetclass(...)` | netclass_exact('x') | exact | pcbnew/pcbexpr_functions.cpp:1961 |
| 59 | function | `A.hasComponentClass(...)` | component_class('x') | exact | pcbnew/pcbexpr_functions.cpp:1962 |
| 60 | function | `A.customProperty(...)` | custom_property('k') (value function, used inside a comparison) | exact | pcbnew/pcbexpr_functions.cpp:1964, :1823-1864 |
| 61 | function | `A.hasCustomProperty(...)` | has_custom_property('k') | exact | pcbnew/pcbexpr_functions.cpp:1965 |

KiCad object classes (inheritance and masked properties are in the KiCad property registrations):

| # | KiCad class | Canonical object type ((+) = not yet in design-graph.md) | Provenance |
|---|---|---|---|
| 62 | `BOARD_CONNECTED_ITEM` | any connected item (net-bearing) | pcbnew/board_connected_item.cpp:218† |
| 63 | `BOARD_ITEM` | any board item | pcbnew/board_item.cpp:538† |
| 64 | `FOOTPRINT` | footprint | pcbnew/footprint.cpp:5260-5261† |
| 65 | `PAD` | pad | pcbnew/pad.cpp:3604-3605† |
| 66 | `PCB_3D_BODY` | body3d (+) | pcbnew/pcb_3d_body.cpp:194† |
| 67 | `PCB_ARC` | track.arc | pcbnew/pcb_track.cpp:3652† |
| 68 | `PCB_BARCODE` | barcode (+) | pcbnew/pcb_barcode.cpp:1015† |
| 69 | `PCB_CONSTRAINT` | none: geometric constraint object (+) | pcbnew/constraints/pcb_constraint.cpp:485† |
| 70 | `PCB_DIMENSION_BASE` | dimension (+) | pcbnew/pcb_dimension.cpp:2122-2124† |
| 71 | `PCB_DIM_ALIGNED` | dimension (+) | pcbnew/pcb_dimension.cpp:2344-2347† |
| 72 | `PCB_DIM_CENTER` | dimension (+) | pcbnew/pcb_dimension.cpp:2501-2504† |
| 73 | `PCB_DIM_LEADER` | dimension.leader (+) | pcbnew/pcb_dimension.cpp:2462-2465† |
| 74 | `PCB_DIM_ORTHOGONAL` | dimension (+) | pcbnew/pcb_dimension.cpp:2387-2391† |
| 75 | `PCB_DIM_RADIAL` | dimension (+) | pcbnew/pcb_dimension.cpp:2419-2422† |
| 76 | `PCB_DRILL_CHART` | table (+) | pcbnew/pcb_drill_chart.cpp:508† |
| 77 | `PCB_DRILL_MAP` | graphic (+) | pcbnew/pcb_drill_map.cpp:391† |
| 78 | `PCB_FIELD` | text.field (+) | pcbnew/pcb_field.cpp:379-381† |
| 79 | `PCB_GENERATED_TABLE` | table (+) | pcbnew/pcb_generated_table.cpp:1112† |
| 80 | `PCB_GENERATOR` | generator (+) | pcbnew/pcb_generator.cpp:214† |
| 81 | `PCB_GRID_ITEM` | graphic (+) | pcbnew/pcb_grid_item.cpp:515† |
| 82 | `PCB_GROUP` | group (+) | pcbnew/pcb_group.cpp:622-623† |
| 83 | `PCB_MARKER` | none: DRC marker, not design data | pcbnew/pcb_marker.cpp:564-565† |
| 84 | `PCB_POINT` | point (+) | pcbnew/pcb_point.cpp:360† |
| 85 | `PCB_REFERENCE_IMAGE` | none: reference image, not design data | pcbnew/pcb_reference_image.cpp:461† |
| 86 | `PCB_SHAPE` | graphic (+) | pcbnew/pcb_shape.cpp:2039-2040† |
| 87 | `PCB_STACKUP_TABLE` | table (+) | pcbnew/pcb_stackup_table.cpp:104† |
| 88 | `PCB_TABLE` | table (+) | pcbnew/pcb_table.cpp:1106-1107† |
| 89 | `PCB_TABLECELL` | table (+) | pcbnew/pcb_tablecell.cpp:372-377† |
| 90 | `PCB_TARGET` | graphic.target (+) | pcbnew/pcb_target.cpp:244† |
| 91 | `PCB_TEXT` | text (+) | pcbnew/pcb_text.cpp:923-924† |
| 92 | `PCB_TEXTBOX` | text.box (+) | pcbnew/pcb_textbox.cpp:1085-1087† |
| 93 | `PCB_TRACK` | track | pcbnew/pcb_track.cpp:3611† |
| 94 | `PCB_TUNING_PATTERN` | generator.tuning (+) | pcbnew/generators/pcb_tuning_pattern.cpp:3322-3323† |
| 95 | `PCB_VIA` | via | pcbnew/pcb_track.cpp:3656† |
| 96 | `PCB_VIA_STACK` | generator.via_stack (+) | pcbnew/generators/pcb_via_stack.cpp:1624-1625† |
| 97 | `PCB_VIA_STITCH` | generator.via_stitch (+) | pcbnew/generators/pcb_via_stitch.cpp:1707-1708† |
| 98 | `ZONE` | zone | pcbnew/zone.cpp:2387† |

All 343 property registrations at this commit map by one rule, `A.<Name>` -> `prop('<Display Name>')`
(underscores become spaces, lookup is case-insensitive, pcbnew/pcbexpr_evaluator.cpp:804-841); the ones with a
first-class canonical predicate are listed in Appendix A.

**Altium query terms** (seen in records or documentation):

| # | Altium term | Canonical scope expression | Quality | Provenance |
|---|---|---|---|---|
| Q1 | `All` | all | exact | golden-rules/odrive-v3.txt:3 |
| Q2 | `InNet('<net>')` | net('n') | exact | W6:1; W9 rec SupplyNets |
| Q3 | `InNetClass('<class>')` | netclass('c') | exact | W5:8; pcbnew/pcb_io/altium/altium_pcb.cpp:127 |
| Q4 | `InDifferentialPair('<pair>')` | diffpair('p') | lossy: Altium pairs are named objects; KiCad infers pairs from _P/_N names | W5:7 |
| Q5 | `InComponent('<ref>')` | component('ref') | exact | golden-rules/odrive-v3.txt:13; W7 rec BoardOutlineClearance |
| Q6 | `InComponentClass('<class>')` | component_class('c') | exact | W8 rec RoomDefinition |
| Q7 | `InPadClass('<class>')` | pad_class('c') | lossy to KiCad: no pad classes (expand to explicit pads or a custom property) | golden-rules/odrive-v3.txt:13 |
| Q8 | `HasFootprint('<name>')` | footprint('name') | lossy: Altium footprint name vs KiCad 'LIB:FP' id | golden-rules/odrive-v3.txt:23, :38 |
| Q9 | `IsPad`, `IsVia`, `IsStitchingVia` | type(pad), type(via), type(via.stitching) | lossy: stitching-via membership has no KiCad property | golden-rules/odrive-v3.txt:12; W6:31 |
| Q10 | `IsBGA`, `IsLCC`, `IsSOIC` | package_style('BGA' \| 'LCC' \| 'SOIC') | lossy: needs a package class from the parts store; KiCad has none | golden-rules/odrive-v3.txt:16-17; W5:26-30 |
| Q11 | `CompPinCount < 5` | pin_count() < 5 | exact to canonical; lossy to KiCad: KiCad registers no pad-count property, so the exporter expands to an explicit refdes list | golden-rules/odrive-v3.txt:18 |
| Q12 | `OnLayer('<layer>')` | on_layer('l') | lossy for user-named and mid layers (for example `OnLayer('Signal 40 1')`, W5:4): mapping a name to `cu:<k>` needs the Board6 layer stack, not decoded (Q-A11); the native layer name is kept in `native.text`. Top and bottom layers map exactly | W5:2-6; W7 rec BoardOutlineClearance |
| Q13 | `OnTopLayer`, `OnBottomLayer`, `OnMultiLayer` | on_layer(top), on_layer(bottom), type(pad.pth) or type(via) | lossy: OnMultiLayer read as through-hole objects (GK) | W6:6 |
| Q14 | `InPolygon`, `InPoly`, `IsPolygon`, `IsPoly` | type(zone) | lossy: Altium polygon primitives vs KiCad zones | pcbnew/pcb_io/altium/altium_parser_pcb.h:1019-1022; pcbnew/pcb_io/altium/altium_parser_pcb.cpp:59-64 |
| Q15 | bare net name, e.g. `DDR3_CKE` | net('n') | lossy: bare-name shorthand unverified (GK) | W5:7 |
| Q16 | `Name = '<pad>'` | on pads: prop('Pad Number') == 'x' | lossy: KiCad's 'Name' property exists only on zones and groups (pcbnew/zone.cpp:2493†, pcbnew/pcb_group.cpp:631†), the pad designator is 'Pad Number' (pcbnew/pad.cpp:3658†); Altium `Name` on other object types has no single KiCad property (GK) | AD13 (eevblog: `IsPad AND (Name = 'Free-0')`) |
| Q17 | `WithinRoom`, `TouchesRoom` | in_area(room, enclosed) / in_area(room, intersects) | exact (no record seen) | AD11 |

Altium query operators (KiCad's tokenizer of Altium queries; nothing in the checkout calls it,
pcbnew/pcb_io/altium/altium_rule_transformer.h:26-67):

| # | Altium operator | Canonical | Quality | Provenance |
|---|---|---|---|---|
| 1 | `And, and, &&` | and | exact | pcbnew/pcb_io/altium/altium_rule_transformer.cpp:45-203 |
| 2 | `Or, or, \|\|` | or | exact | pcbnew/pcb_io/altium/altium_rule_transformer.cpp:45-203 |
| 3 | `Xor` | (a and not b) or (not a and b) | exact | pcbnew/pcb_io/altium/altium_rule_transformer.cpp:45-203 |
| 4 | `Not` | not | exact | pcbnew/pcb_io/altium/altium_rule_transformer.cpp:45-203 |
| 5 | `=, <>` | ==, != | exact | pcbnew/pcb_io/altium/altium_rule_transformer.cpp:45-203 |
| 6 | `<, <=, >, >=` | same | exact | pcbnew/pcb_io/altium/altium_rule_transformer.cpp:45-203 |
| 7 | `+, -, *, /, Div, Mod` | arithmetic (Div = integer division, Mod = remainder (GK)) | lossy: Div and Mod semantics unverified (GK) | pcbnew/pcb_io/altium/altium_rule_transformer.cpp:45-203 |
| 8 | `Between` | x >= lo and x <= hi | lossy: inclusive bounds unverified (GK) | pcbnew/pcb_io/altium/altium_rule_transformer.cpp:45-203 |
| 9 | `Like` | ~= (glob) | lossy: wildcard set unverified (GK) | pcbnew/pcb_io/altium/altium_rule_transformer.cpp:45-203 |
| 10 | `True, False` | true, false | exact | pcbnew/pcb_io/altium/altium_rule_transformer.cpp:45-203 |
| 11 | `'...' strings` | '...' (Altium has no escape; a canonical string containing `\'` is not exportable to Altium, section 2.4) | exact on import | pcbnew/pcb_io/altium/altium_rule_transformer.cpp:45-203, :118 |
| 12 | `integer constants` | count | exact | pcbnew/pcb_io/altium/altium_rule_transformer.cpp:45-203 |
| 13 | `precedence` | And binds tighter than Or (GK); && and \|\| are the lower-precedence forms | lossy: unverified (GK); importer parenthesises | pcbnew/pcb_io/altium/altium_rule_transformer.cpp:45-203 |

Altium legacy scope list (`SCOPE<i>_<j>_KIND`, stored beside the query in ASCII PcbDocs) and object ids:

| # | TScopeKind (low to high precedence) | Canonical scope expression | Provenance |
|---|---|---|---|
| 1 | `Board` | all | unverified (L1: confirm via Altium scripting, Q-A12) |
| 2 | `LayerClass` | on_layer(layer class members) | unverified (L1: confirm via Altium scripting, Q-A12) |
| 3 | `Layer` | on_layer('l') | unverified (L1: confirm via Altium scripting, Q-A12) |
| 4 | `ObjectKind` | type(...) | unverified (L1: confirm via Altium scripting, Q-A12) |
| 5 | `Footprint` | footprint('name') | unverified (L1: confirm via Altium scripting, Q-A12) |
| 6 | `ComponentClass` | component_class('c') | unverified (L1: confirm via Altium scripting, Q-A12) |
| 7 | `Component` | component('ref') | unverified (L1: confirm via Altium scripting, Q-A12) |
| 8 | `NetClass` | netclass('c') | unverified (L1: confirm via Altium scripting, Q-A12) |
| 9 | `Net` | net('n') | unverified (L1: confirm via Altium scripting, Q-A12) |
| 10 | `FromToClass` | from_to over a class (GK) | unverified (L1: confirm via Altium scripting, Q-A12) |
| 11 | `FromTo` | from_to('a', 'b') | unverified (L1: confirm via Altium scripting, Q-A12) |
| 12 | `PadClass` | pad_class('c') | unverified (L1: confirm via Altium scripting, Q-A12) |
| 13 | `PadSpec` | prop comparisons on pad size and shape (GK) | unverified (L1: confirm via Altium scripting, Q-A12) |
| 14 | `ViaSpec` | prop comparisons on via size (GK) | unverified (L1: confirm via Altium scripting, Q-A12) |
| 15 | `FootprintPad` | in_footprint(f) and prop('Pad Number') == p | unverified (L1: confirm via Altium scripting, Q-A12) |
| 16 | `Pad` | component(ref) and prop('Pad Number') == p | unverified (L1: confirm via Altium scripting, Q-A12) |
| 17 | `Region` | in_area(region, intersects) | unverified (L1: confirm via Altium scripting, Q-A12) |

| # | TScopeObjectId | Canonical | Provenance |
|---|---|---|---|
| 1 | `Wire` | type(track) | unverified (L1: confirm via Altium scripting, Q-A12) |
| 2 | `Pin` | type(pad) | unverified (L1: confirm via Altium scripting, Q-A12) |
| 3 | `Smd` | type(pad.smd) | unverified (L1: confirm via Altium scripting, Q-A12) |
| 4 | `Via` | type(via) | unverified (L1: confirm via Altium scripting, Q-A12) |
| 5 | `Fill` | type(graphic.fill) (GK) | unverified (L1: confirm via Altium scripting, Q-A12) |
| 6 | `Polygon` | type(zone) | unverified (L1: confirm via Altium scripting, Q-A12) |
| 7 | `KeepOut` | type(keepout) | unverified (L1: confirm via Altium scripting, Q-A12) |

**Eagle.** Eagle has no query language. Its scopes are: the whole board (every `<param>`), an object-type pair
(`md*` rows), a layer group (`rv*`, `rl*` rows), a net class (`width`, `drill`) and a class pair (the clearance
matrix):

| Eagle scope | Canonical scope | Quality | Provenance |
|---|---|---|---|
| board-wide `<param>` | `expr: all` | exact | pcbnew/pcb_io/eagle/pcb_io_eagle.cpp:162-229 |
| object-type pair (`mdWirePad` ...) | `expr: type(track)`, `expr_b: type(pad.pth)`, `pair: [different_net]` | lossy: Eagle's "pad" vs "smd" split must come from the pad attribute; "wire" includes polygon outlines (GK-E) | OMX; pcbnew/pcb_io/eagle/pcb_io_eagle.cpp:175-224 |
| layer group (pad top / inner / bottom; via outer / inner) | `layers: [cu:1]`, `[inner]`, `[cu:N]`, `[outer]` | exact | pcbnew/pcb_io/eagle/pcb_io_eagle.h:115-127 |
| class `N` (`<signal class="N">`) | `expr: netclass_exact('<name of N>')` | exact | pcbnew/pcb_io/eagle/pcb_io_eagle.cpp:2779-2789 |
| class matrix entry N:M | `expr: netclass_exact(N)`, `expr_b: netclass_exact(M)`, symmetric | exact (KiCad's importer writes `A.NetClass == ... && B.NetClass == ...`, pcbnew/pcb_io/eagle/pcb_io_eagle.cpp:2707-2709) | pcbnew/drc/drc_rule_condition.cpp:105-117 |
| restrict shapes (layers 41-43) | `in_area(<shape>, intersects)` with `disallow` | exact | pcbnew/pcb_io/eagle/pcb_io_eagle.cpp:121-147 |

### 2.6 Lossy cases

1. **KiCad operator precedence.** `||` appears to bind tighter than `&&` (inferred, not run, Q-K6), and `!` binds
   tighter than the comparison operators (grammar.lemon:31-32). Import keeps KiCad's parse tree; export parenthesises
   every mixed and/or and writes `!( ... )` around every negated non-atomic term.
2. **Altium objects without a KiCad property**: pad classes (`InPadClass`), package styles (`IsBGA`, `IsLCC`, `IsSOIC`),
   pin counts (`CompPinCount`), stitching vias (`IsStitchingVia`), named differential pairs. Canonical keeps them
   (`pad_class()`, `package_style()`, the value function `pin_count()`); the KiCad exporter expands them to explicit
   lists (refdes, pads, nets) at export time, or writes a custom property and a value comparison on it
   (`A.customProperty('pad_class') == 'All Pads'`, canonical `custom_property('pad_class') == 'All Pads'`, master
   only; section 2.4 value functions) and marks the rule lossy.
3. **Altium footprint names vs KiCad library ids**: `HasFootprint('1704936')` has no `LIB:` part; the KiCad exporter
   matches on the footprint name with a wildcard library (`'*:1704936'`), which may over-match.
4. **Altium binary rules with SCOPE2 not `All`** export to KiCad as one condition `(<A part>) && (<B part>)`, which is
   exact because KiCad evaluates both orders. The reverse (a KiCad condition mixing A and B in one sub-expression,
   for example `A.NetClass == B.NetClass`) has no separable Altium form; the Altium exporter writes the closest
   separable rules and marks the rule lossy.
5. **Eagle object-type pairs** need the SMD / through-hole split of pads; Eagle's "wire" class is assumed to mean
   signal tracks (GK-E).
6. **Layer names**: Altium mid-layer ordinals need the PcbDoc layer stack (Q-A11); Eagle inner layers need
   `layerSetup` (GK-E).
7. **Same-net clearance**: Altium `SameNetOnly` clearance has no KiCad `clearance` form (KiCad clearance is between
   different nets); it maps to `physical_clearance` with `pair: [same_net]`, which KiCad writes as
   `physical_clearance` with `A.NetName == B.NetName`.
8. **Item overrides** (pad clearance, mask, paste, zone connection) are object data in KiCad and Eagle but rules in
   the canonical model (origin `item_override`). Altium "Rule" vs "Manual" mask modes on pads (KiCad's importer
   treats Rule as inherit, pcbnew/pcb_io/altium/altium_pcb.cpp:5224-5228) round-trip only if the mode is kept in
   `params`.

### 2.7 Severity

| Canonical | KiCad | Altium | Eagle |
|---|---|---|---|
| `blocker` | none (kernel only: release-blocking) | none | none |
| `error` | `(severity error)`; severity map value `error` | DRC violation (GK); ERC levels Error and Fatal Error (AD14) | DRC error (GK-E) |
| `warning` | `(severity warning)`; `warning` | ERC level Warning (AD14) | none (GK-E) |
| `info` | none | none | none |
| `ignore` | `(severity ignore)`: matched objects are not reported (pcbnew/drc/drc_test_provider_copper_clearance.cpp:281); severity map value `ignore` | rule disabled (`ENABLED=FALSE`); ERC level No Report (AD14) | check switched off (`check*` = 0, GK-E) |
| `inherit` | no `(severity ...)` clause: the board-level severity of the violation code applies (pcbnew/drc/drc_rule.cpp:44; pcbnew/pcb_marker.cpp:427-430†) | every PCB rule (no per-rule severity field in any record) | every rule |
| tag `waive` | `(severity exclusion)`: markers counted as excluded (pcbnew/pcb_marker.cpp:427-428†) | none | none |
| finding status `dismissed` | `drc_exclusions`, `erc_exclusions` | DRC/ERC waivers (GK) | approved errors (GK-E) |

`inherit` resolves through the rule set's per-kind default severity, which importers fill from KiCad's
`rule_severities` maps (section 4.1.12, 4.1.13) and Altium's `[Electrical Rules Check]` levels (once the type
numbers are mapped, Q-A6). A KiCad violation code with a fixed severity (`padstack_invalid`, `generic_error`,
pin-to-pin matrix codes) keeps it (pcbnew/drc/drc_item.cpp:486-491; eeschema/erc/erc_settings.cpp:322-341).

### 2.8 Fit with `schemas/rule.schema.json`

The schema before this work could not express a geometric rule: it had no kind, no object-pair scope, no layer
filter and no min / opt / max. It was extended in the same change, minimally and backward-compatibly: new optional
fields `kind`, `variant`, `scope.expr`, `scope.expr_b`, `scope.pair`, `scope.symmetric`, `layers`, `constraint`
(`min`, `opt`, `max`, `unit`, `value`, `options`), `priority`, `origin`, `floor_exempt`, `enabled`, `native`; and two
new severity values, `ignore` and `inherit`. Task F1 added the `origin` value `item_default` (section 5, P1 and P9). All
previously required fields stay required; no existing property
changed type. A rule record from an import therefore carries `check` (`pcbbrain_checks.drc.<kind>`), `title`,
`provenance` (the source file and the native rule name) and `confidence` (`high` for parsed tool rules) as before.

## 3. Canonical taxonomy

129 kinds in fourteen categories. Counts of round-trip quality per tool:

| Category | Kinds | KiCad exact / lossy / none | Altium exact / lossy / none | Eagle exact / lossy / none |
|---|---|---|---|---|
| Electrical and connectivity | 11 | 2 / 6 / 3 | 4 / 5 / 2 | 0 / 1 / 10 |
| Clearance | 9 | 5 / 1 / 3 | 3 / 5 / 1 | 1 / 1 / 7 |
| Geometry, width and size | 5 | 4 / 0 / 1 | 1 / 1 / 3 | 2 / 0 / 3 |
| Holes and vias | 11 | 5 / 5 / 1 | 5 / 2 / 4 | 3 / 4 / 4 |
| Solder mask, paste, silkscreen and text | 12 | 3 / 8 / 1 | 5 / 2 / 5 | 3 / 2 / 7 |
| Copper pour, plane and thermal relief | 8 | 6 / 2 / 0 | 4 / 1 / 3 | 1 / 3 / 4 |
| Routing style and topology | 13 | 2 / 2 / 9 | 6 / 4 / 3 | 0 / 1 / 12 |
| High speed: length, differential pair, skew | 13 | 8 / 3 / 2 | 4 / 6 / 3 | 0 / 2 / 11 |
| Placement, courtyard, height and keepout | 9 | 2 / 6 / 1 | 2 / 4 / 3 | 0 / 1 / 8 |
| Testpoint and assembly | 3 | 0 / 0 / 3 | 0 / 2 / 1 | 0 / 1 / 2 |
| Manufacturing (fab capability) | 5 | 1 / 2 / 2 | 0 / 1 / 4 | 0 / 1 / 4 |
| Signal integrity | 9 | 0 / 0 / 9 | 1 / 8 / 0 | 0 / 0 / 9 |
| Schematic and ERC | 14 | 14 / 0 / 0 | 0 / 14 / 0 | 0 / 1 / 13 |
| General: custom assertions, data integrity, engine settings (added category) | 7 | 2 / 4 / 1 | 1 / 1 / 5 | 1 / 1 / 5 |
| **Total** | 129 | 54 / 39 / 36 | 36 / 56 / 37 | 11 / 19 / 99 |

Reading the counts: KiCad covers the electrical, clearance, geometry, hole, high-speed and ERC kinds; 18 of its 39
"lossy" entries are built-in checks with no constraint of their own: 16 that KiCad switches on or off only board-wide,
and 2 (short_circuit, tracks_crossing) that can be scoped only indirectly through the pair's clearance rule. Altium is the only
source for routing directives, SMT entry rules, testpoints and signal integrity; 45 of its 56 "lossy" entries are
kinds whose parameter keys, meanings or ERC type numbers no source gave (they become exact once L1 confirms them).
Eagle maps to 30 kinds, 11 exactly; 13 of its 19 losses are parameter meanings that no official source confirms.

### 3.1 Electrical and connectivity

| Canonical id | Meaning | Parameters (units) | Arity | KiCad | Altium | Eagle | Mapping quality | Provenance |
|---|---|---|---|---|---|---|---|---|
| `short_circuit` | Copper of two different nets touches or overlaps. | allowed (bool, default false) | B | DRC code `shorting_items`; no constraint of its own: it is tested only when the pair's resolved `clearance` is > 0 and its severity is not `ignore` (item, pad, zone and zone-to-zone paths; knockout text against a zone is tested at >= 0), so a clearance rule with min <= 0 or `(severity ignore)` turns shorts off for its scope | ShortCircuit (`ALLOWED`) | none | KiCad: lossy: scoped only indirectly through the pair's clearance rule; no `allowed` value (import: a KiCad clearance rule with min <= 0 mm or severity ignore implies short_circuit.allowed=true for its scope; export of a scoped allowed=true: a clearance rule with min 0 mm on that scope, which also drops the clearance check there). Altium: exact. Eagle: none | pcbnew/drc/drc_item.h:49; pcbnew/drc/drc_test_provider_copper_clearance.cpp:275-281, :315, :495, :531, :601, :866-879, :1311; golden-rules/odrive-v3.txt:37 |
| `tracks_crossing` | Two tracks of different nets cross on one layer. | none | B | DRC code `tracks_crossing`; tested inside the same gate as shorts (pair clearance > 0, severity not `ignore`) | none (part of ShortCircuit) | none | KiCad: lossy: scoped only indirectly through the pair's clearance rule, as short_circuit. Altium: lossy: covered by ShortCircuit, not separable. Eagle: none | pcbnew/drc/drc_item.h:54; pcbnew/drc/drc_test_provider_copper_clearance.cpp:281-294 |
| `unrouted_net` | A net still has unrouted connections (open ratsnest). | check_bad_connections (bool) | U (net) | DRC code `unconnected_items` | UnRoutedNet (`CHECKBADCONNECTIONS`) | none | KiCad: lossy: board-wide severity only, not scopable. Altium: exact. Eagle: none | pcbnew/drc/drc_item.h:48; golden-rules/odrive-v3.txt:36; W6:27 |
| `unconnected_pin` | A pad on the board has no connection at all. | none | U (pad) | none on the board (schematic side: `erc_unconnected_pin`) | UnConnectedPin (keys unknown) | none | KiCad: none. Altium: lossy: keys and meaning unverified (GK). Eagle: none | AD4; unverified (GK) |
| `dangling_track` | A track end connects to nothing (antenna). | max_stub (mm; default 0) | U (track) | DRC code `track_dangling` | NetAntennae (`NETANTENNAETOLERANCE`) | `checkWireStubs` (flag) | KiCad: lossy: board-wide severity only, not scopable; tolerance fixed at 0. Altium: exact. Eagle: lossy: on/off flag only; meaning unverified (GK) | pcbnew/drc/drc_item.h:60; golden-rules/odrive-v3.txt:9; golden-rules/moteus.txt:91 |
| `dangling_via` | A via connects to nothing, or to copper on one layer only. | none | U (via) | `via_dangling` (severity-only constraint); DRC code `via_dangling` | none found | none | KiCad: exact. Altium: none. Eagle: none | common/drc_rules.keywords:73; pcbnew/drc/drc_test_provider_connectivity.cpp:94-97 |
| `isolated_copper` | A copper island is not connected to the rest of its net. | min_area (mm2) | U | DRC code `isolated_copper` | DRC option `REPORTDEADCOPPER` with `DEADCOPPERMINAREA` | none | KiCad: lossy: board-wide severity only, not scopable; no area threshold. Altium: lossy: DRC option, not a rule; area unit unresolved. Eagle: none | pcbnew/drc/drc_item.h:57; odrive-DRC D5 |
| `unpoured_zone` | A copper pour is unfilled, or was modified after its last fill. | none | U (zone) | none (no DRC code at this commit) | UnpouredPolygon (UI: Modified Polygon) | none | KiCad: none. Altium: exact (meaning unverified (GK)). Eagle: none | golden-rules/odrive-v3.txt:8; W13:95 |
| `plane_integrity` | Plane layers are unbroken and holes sit correctly against plane voids. | variant: broken_plane, internal_plane_warning, straddling_hole, hole_in_void | U | none | DRC options `REPORTBROKENPLANES`, `INTERNALPLANEWARNINGS`, `REPORTSTRADLINGHOLES`, `REPORTHOLESINVOIDS` | none | KiCad: none. Altium: lossy: DRC options, meanings unverified (GK). Eagle: none | odrive-DRC D4, D7 |
| `net_tie` | Nets intentionally joined inside a net-tie footprint; clearance and short checks between them are waived there. | pad_groups (list of pad-number sets) | B | footprint `net_tie_pad_groups` (item data); clearance forced to 0 mm | none found (unverified (GK)) | none | KiCad: exact. Altium: none. Eagle: none | pcbnew/drc/drc_engine.cpp:1142-1176 |
| `board_schematic_parity` | Board and schematic disagree on components, nets, attributes or fields. | variant: missing_footprint, extra_footprint, duplicate_footprint, net_conflict, symbol_mismatch, footprint_filter_mismatch, field_mismatch | U | 7 DRC codes of the `schematic_parity` provider | none as rules (`[Difference Levels]`, meaning unverified (GK)) | none | KiCad: lossy: board-wide severity only, not scopable. Altium: lossy: comparator settings, meanings unverified (GK). Eagle: none | pcbnew/drc/drc_test_provider_schematic_parity.cpp:31-40; golden/odrive-v3/src/v3/two_ax_PCB.PrjPcb:4378 |

### 3.2 Clearance

| Canonical id | Meaning | Parameters (units) | Arity | KiCad | Altium | Eagle | Mapping quality | Provenance |
|---|---|---|---|---|---|---|---|---|
| `clearance` | Minimum distance between copper objects of different nets on one layer (electrical clearance). | min (mm; in KiCad a value <= 0 mm disables the clearance, short and crossing tests for the pair, section 3.1) | B | `clearance`; floor `rules.min_clearance`; netclass `clearance`; item local clearance | Clearance (`GAP`, `NETSCOPE`, `OBJECTCLEARANCES`, `IGNOREPADTOPADCLEARANCEINFOOTPRINT`) | `md*` object-pair distances; class clearance matrix | KiCad: exact. Altium: exact (each `OBJECTCLEARANCES` entry becomes its own rule). Eagle: lossy: object-type pairs only; SMD vs through-hole pad needs the pad attribute; Eagle's combination of matrix and md* unverified (GK) | common/drc_rules.keywords:8; pcbnew/drc/drc_engine.cpp:2086-2162; golden-rules/odrive-v3.txt:35; golden-rules/moteus.txt:91, :93 |
| `physical_clearance` | Minimum distance between any two objects on a layer regardless of net, non-copper layers included. | min (mm) | B | `physical_clearance` (deprecated `mechanical_clearance`); implicit barcode rule | Clearance with `NETSCOPE` `AnyNet` or `SameNetOnly` | none | KiCad: exact. Altium: lossy: copper objects only. Eagle: none | common/drc_rules.keywords:48; pcbnew/drc/drc_engine.cpp:253-260; W6:1 |
| `hole_clearance` | Minimum distance from a drilled hole edge to copper of another net. | min (mm) | B | `hole_clearance`; floor `rules.min_hole_clearance` | Clearance `OBJECTCLEARANCES` entries naming `Hole` | none (see Q-E3 on `mdDrill`) | KiCad: exact. Altium: lossy: only as an object-pair entry; Hole semantics unverified (GK). Eagle: none | common/drc_rules.keywords:24; pcbnew/drc/drc_engine.cpp:223-226; W5:1 |
| `physical_hole_clearance` | Minimum distance from a round drilled hole to any object regardless of net (drilling tolerance). | min (mm) | B | `physical_hole_clearance` (deprecated `mechanical_hole_clearance`) | none | none | KiCad: exact. Altium: none. Eagle: none | common/drc_rules.keywords:49 |
| `edge_clearance` | Minimum distance from copper (and milled pad holes) to the board edge. | min (mm) | B (object, edge) | `edge_clearance`; floor `rules.min_copper_edge_clearance` | BoardOutlineClearance (`GAP`, `OBJECTCLEARANCES`) | `mdCopperDimension` | KiCad: exact. Altium: exact. Eagle: exact (meaning from OMX, third-party) | common/drc_rules.keywords:17; W6:2; golden-rules/moteus.txt:91 |
| `creepage` | Minimum surface (creepage) path between conductors of two nets. | min (mm); min_groove_width (mm); ignore_internal_layers (bool) | B | `creepage`; `rules.min_groove_width` | Creepage (keys unknown) | none | KiCad: lossy: no ignore_internal_layers option. Altium: lossy: keys unknown (AD5 summary only). Eagle: none | common/drc_rules.keywords:13; pcbnew/drc/drc_test_provider_creepage.cpp:308; AD5 |
| `z_axis_clearance` | Minimum distance between copper on different layers, through the dielectric. | min (mm) | B | none (KiCad clearance is same-layer) | ZAxisClearance (keys unknown) | none | KiCad: none. Altium: lossy: keys unknown. Eagle: none | AD12 |
| `plane_clearance` | Clearance from a power plane to holes and pads of other nets. | min (mm) | B | none (export approximation: zone local clearance on plane-layer zones) | PlaneClearance (`CLEARANCE`) | none | KiCad: none. Altium: exact. Eagle: none | golden-rules/odrive-v3.txt:33; pcbnew/pcb_io/altium/altium_pcb.cpp:3612-3636 |
| `wire_bonding` | Altium AD26 wire-bonding rule (bond-wire clearance per AD12); keys and exact meaning unknown. | unknown | unknown | none | WireBonding (no record) | none | KiCad: none. Altium: lossy: meaning and keys unknown. Eagle: none | AD12 |

### 3.3 Geometry, width and size

| Canonical id | Meaning | Parameters (units) | Arity | KiCad | Altium | Eagle | Mapping quality | Provenance |
|---|---|---|---|---|---|---|---|---|
| `track_width` | Width of track segments and arcs. | min / opt / max (mm); opt = router default | U (track) | `track_width`; floor `rules.min_track_width`; netclass `track_width`, `diff_pair_width` (opt); tuning profile | Width (`MINLIMIT`, `PREFEREDWIDTH`, `MAXLIMIT`, per-layer keys); DiffPairsRouting per-layer widths | `msWidth` (min); class `width` (min per class) | KiCad: exact. Altium: exact (one rule per layer; mid-layer ordinals need the layer stack). Eagle: exact (min only) | common/drc_rules.keywords:68; golden-rules/odrive-v3.txt:31; golden-rules/moteus.txt:91 |
| `connection_width` | Minimum width of copper along a connection, including necks inside pours. | min (mm) | B (A, and optionally B: the first two items that contribute copper at the narrow point, e.g. pad and zone; section 2.2) | `connection_width`; floor `rules.min_connection` | none found (Width `CHECKCONNECTEDCOPPER`, meaning unverified (GK)) | none | KiCad: exact. Altium: none. Eagle: none | common/drc_rules.keywords:10; pcbnew/drc/drc_engine.cpp:177-179; pcbnew/drc/drc_test_provider_connection_width.cpp:450-454; pcbnew/drc/drc_rule.cpp:88-118 (not unary) |
| `track_segment_length` | Length of each track segment or arc. | min / max (mm) | U (track) | `track_segment_length` | none | none | KiCad: exact. Altium: none. Eagle: none | common/drc_rules.keywords:69 |
| `track_angle` | Angle between two connected straight track segments. | min / max (deg) | B (segment pair) | `track_angle` | AcuteAngle (min; keys unknown) | none | KiCad: exact. Altium: lossy: min only; keys unknown. Eagle: none | common/drc_rules.keywords:67; AD2 |
| `pad_shape_rule` | Generation: pad shape when the library gives none (shape overrides, elongation, SMD corner rounding). | shape_top / shape_bottom / shape_first (enum); elongation_long / elongation_offset (%); roundness (ratio); roundness_min / roundness_max (mm) | U (pad) | none (geometry is baked into footprints) | none | `psTop`, `psBottom`, `psFirst`, `psElongationLong`, `psElongationOffset`, `srRoundness`, `srMinRoundness`, `srMaxRoundness` | KiCad: none. Altium: none. Eagle: exact (shape codes unverified (GK)) | pcbnew/pcb_io/eagle/pcb_io_eagle.cpp:191-202; golden-rules/moteus.txt:91 |

### 3.4 Holes and vias

| Canonical id | Meaning | Parameters (units) | Arity | KiCad | Altium | Eagle | Mapping quality | Provenance |
|---|---|---|---|---|---|---|---|---|
| `hole_size` | Drill diameter of pad and via holes (oval holes: minor axis against min, major against max). | min / opt / max (mm); or min_pct / max_pct (% of pad size) | U (hole) | `hole_size` (deprecated `hole`); floors `rules.min_through_hole_diameter`, `rules.min_microvia_drill`; netclass `via_drill`, `microvia_drill` (opt) | HoleSize (`ABSOLUTEVALUES`, `MINLIMIT`, `MAXLIMIT`, `MINPERCENT`, `MAXPERCENT`); RoutingVias hole keys | `msDrill`, `msMicroVia` (min); class `drill` (min per class) | KiCad: lossy: no percent-of-pad form. Altium: exact. Eagle: exact (min only) | common/drc_rules.keywords:25; golden-rules/odrive-v3.txt:20, :27; golden-rules/moteus.txt:91 |
| `via_diameter` | Copper diameter of vias. | min / opt / max (mm) | U (via) | `via_diameter`; floors `rules.min_via_diameter`, `rules.min_microvia_diameter`; netclass `via_diameter`, `microvia_diameter` (opt) | RoutingVias (`WIDTH`, `MINWIDTH`, `MAXWIDTH`) | none (derived as drill + 2 x ring, see annular_ratio) | KiCad: exact. Altium: exact. Eagle: none | common/drc_rules.keywords:74; golden-rules/odrive-v3.txt:27 |
| `annular_width` | Width of the copper ring around a plated hole. | min / max (mm) | U (pad, via) | `annular_width`; floor `rules.min_via_annular_width` | MinimumAnnularRing (`MINIMUMRING`) | `rlMin*`, `rlMax*` per layer group (pad top / inner / bottom, via outer / inner, micro-via outer / inner) | KiCad: exact. Altium: exact (min only). Eagle: exact (layer groups become layer filters) | common/drc_rules.keywords:1; golden-rules/odrive-v3.txt:3; golden-rules/moteus.txt:91 |
| `annular_ratio` | Generation: ring = clamp(min, ratio x drill, max) for pads and vias without an explicit diameter. | ratio (fraction) per layer group; min / max (mm); honour_library_diameter (flags) | U (pad, via) | none | none | `rv*` (7 keys); `useDiameter` | KiCad: none. Altium: none. Eagle: exact (`useDiameter` bit layout unverified (GK)) | pcbnew/pcb_io/eagle/pcb_io_eagle.cpp:2084-2091; golden-rules/moteus.txt:91 |
| `hole_to_hole` | Minimum distance between the edges of two mechanically drilled holes. | min (mm) | B | `hole_to_hole`; floor `rules.min_hole_to_hole` | HoleToHoleClearance (`GAP`) | `mdDrill` (per OMX; see Q-E3) | KiCad: exact. Altium: exact. Eagle: lossy: meaning from a third-party glossary only | common/drc_rules.keywords:26; golden-rules/odrive-v3.txt:14; golden-rules/moteus.txt:91 |
| `hole_colocated` | Two holes at the same position. | none | B | DRC code `holes_co_located` | none | none | KiCad: lossy: board-wide severity only, not scopable. Altium: none. Eagle: none | pcbnew/drc/drc_item.h:62 |
| `hole_aspect_ratio` | Maximum ratio of hole depth to drill diameter. | max (ratio) | U (hole) | `hole_aspect_ratio` (master only) | none | `msBlindViaRatio` (min drill / depth, blind vias) | KiCad: exact. Altium: none. Eagle: lossy: inverted ratio, blind vias only | common/drc_rules.keywords:23; golden-rules/moteus.txt:91 |
| `microvia_stack` | Stacked micro-via limits and stack integrity. | variant: depth, malformed_span, not_filled, crosses_core; max_depth (count); allow_stacked (bool) | U (via) | `microvia_stack_depth`; DRC codes `malformed_microvia_stack_span`, `microvia_stack_not_filled`, `microvia_crosses_core` | HoleToHoleClearance `ALLOWSTACKEDMICROVIAS` | none | KiCad: exact. Altium: lossy: allow or forbid only. Eagle: none | common/drc_rules.keywords:38; pcbnew/drc/drc_test_provider_via_stack.cpp:79; golden-rules/odrive-v3.txt:14 |
| `via_span` | Which via types and drill layer pairs are allowed. | allowed_types (set: through, blind, buried, micro); allowed_pairs (list of layer pairs); enforce_pairs (bool) | U (via) | `disallow` kinds `through_via`, `blind_via`, `buried_via`, `micro_via` | LayerPairs (`ENFORCE`); RoutingVias `VIASTYLE` | `layerSetup` (drillable spans) | KiCad: lossy: no drill-pair list. Altium: lossy: `VIASTYLE` values other than `Through Hole` unverified (GK). Eagle: lossy: `layerSetup` syntax unverified (GK) | pcbnew/drc/drc_rule_parser.cpp:630-671; golden-rules/odrive-v3.txt:21, :27; golden-rules/moteus.txt:91 |
| `vias_under_smd` | Whether vias may sit inside SMD pads. | allowed (bool) | U (via) | `disallow via` with a pad-intersection condition (graphical editor output) | ViasUnderSMD (`ALLOWED`) | `mnLayersViaInSmd` (meaning unknown) | KiCad: lossy: needs a geometric condition. Altium: exact. Eagle: lossy: meaning unverified (GK) | pcbnew/drc/rule_editor/drc_re_vias_under_smd_constraint_data.h:86-101; W7 rec ViasUnderSMD |
| `post_machining_conflict` | A track sits on a copper layer removed by back-drilling or post-machining. | none | U (track) | DRC code `track_on_post_machined_layer` | none | none | KiCad: lossy: board-wide severity only, not scopable. Altium: none. Eagle: none | pcbnew/drc/drc_item.h:133 |

### 3.5 Solder mask, paste, silkscreen and text

| Canonical id | Meaning | Parameters (units) | Arity | KiCad | Altium | Eagle | Mapping quality | Provenance |
|---|---|---|---|---|---|---|---|---|
| `solder_mask_expansion` | Growth of the solder-mask aperture around pads, vias, tracks and shapes. | opt (mm, may be negative); or ratio (fraction of the smaller pad side) with min / max clamp (mm) | U | `solder_mask_expansion`; board `pad_to_mask_clearance`; item `solder_mask_margin` | SolderMaskExpansion (`EXPANSION`) | `mvStopFrame` with `mlMinStopFrame`, `mlMaxStopFrame` | KiCad: lossy: no ratio form (exact when absolute, or min = max). Altium: exact. Eagle: exact | common/drc_rules.keywords:54; golden-rules/odrive-v3.txt:25; golden-rules/moteus.txt:91 |
| `solder_mask_sliver` | Minimum width of the solder-mask web between two apertures. | min (mm) | B | board `solder_mask_min_width` (used only as the bridge criterion, see solder_mask_bridge); `solder_mask_sliver` (parsed, no consumer found) | MinimumSolderMaskSliver (`MINSOLDERMASKWIDTH`) | none | KiCad: lossy: KiCad never reports a sliver. Its board value only inflates apertures by half the web, so different-net apertures closer than the web are reported as `solder_mask_bridge`, while a narrow web between same-net apertures is not reported at all. Altium: exact. Eagle: none | pcbnew/drc/drc_test_provider_solder_mask.cpp:205, :506, :584-587, :1003; golden-rules/odrive-v3.txt:13 |
| `solder_mask_bridge` | Mask apertures of different nets merge into one opening. | allowed (bool); allow_in_footprint (bool); min_web (mm: apertures closer than this count as merged); min_copper_to_aperture (mm); aperture_alone (bool option: a bare mask aperture matching the scope is exempt on its own) | B | DRC code `solder_mask_bridge`; `bridged_mask` (unary, severity-only exemption); board `solder_mask_min_width` (the web criterion; the provider reports bridges against an internal rule named "board setup solder mask min width"); `allow_soldermask_bridges_in_footprints`; `rules.solder_mask_to_copper_clearance` | none found | `checkStop` (flag; meaning unverified (GK)) | KiCad: lossy: `bridged_mask` is unary and exempts a copper-to-copper bridge only when both participants match with severity `ignore` (a bare aperture alone suffices); import writes expr = expr_b = condition with aperture_alone=true (section 4.1.1 row 26). Altium: none. Eagle: lossy: flag only, meaning unverified (GK) | common/drc_rules.keywords:5; pcbnew/drc/drc_rule.cpp:111; pcbnew/drc/drc_test_provider_solder_mask.cpp:77, :517, :584-610, :840, :1003 |
| `via_protection` | Default via tenting, covering, plugging, capping and filling. | tent_top / tent_bottom, cover_top / cover_bottom, plug_top / plug_bottom, cap, fill (bool); tent_max_drill (mm) | U (via) | board `tenting`, `covering`, `plugging`, `capping`, `filling`; per-via properties | SolderMaskExpansion `ISTENTINGTOP`; mask rules scoped `IsVia` | `mlViaStopLimit`; via `alwaysstop` | KiCad: lossy: no drill threshold. Altium: lossy: tenting is expressed through mask-expansion rules. Eagle: exact | include/board_design_settings.h:701-714; W5:16; W6:31; golden-rules/moteus.txt:91 |
| `solder_paste_margin` | Change of the paste aperture relative to the pad. | abs (mm, usually negative); ratio (fraction); clamp_min / clamp_max (mm) | U (pad) | `solder_paste_abs_margin`, `solder_paste_rel_margin`; board `pad_to_paste_clearance`, `pad_to_paste_clearance_ratio`; item margins | PasteMaskExpansion (`EXPANSION`) | `mvCreamFrame` with `mlMinCreamFrame`, `mlMaxCreamFrame` | KiCad: lossy: no clamps; DRU encoding of the ratio unresolved (Q-K1). Altium: exact (absolute only). Eagle: exact | common/drc_rules.keywords:55-56; golden-rules/odrive-v3.txt:26; golden-rules/hackrf-one.txt:16 |
| `silk_clearance` | Minimum distance between silkscreen objects. | min (mm) | B | `silk_clearance` (key `silk_overlap`); floor `rules.min_silk_clearance` | SilkToSilkClearance (`SILKTOSILKCLEARANCE`) | none | KiCad: exact. Altium: exact. Eagle: none | common/drc_rules.keywords:52; golden-rules/odrive-v3.txt:11 |
| `silk_to_mask_clearance` | Minimum distance from silkscreen to a mask aperture or to exposed copper. | min (mm); measure_to (aperture, exposed_copper) | B | `silk_clearance` checked against mask items (key `silk_over_copper`) | SilkToSolderMaskClearance (`MINSILKSCREENTOMASKGAP`, `CLEARANCETOEXPOSEDCOPPER`) | none | KiCad: lossy: shares the silk_clearance value unless split by a layer condition; no measure_to. Altium: exact. Eagle: none | pcbnew/drc/drc_test_provider_solder_mask.cpp:329; golden-rules/odrive-v3.txt:12 |
| `silk_to_edge_clearance` | Minimum distance from silkscreen to the board edge or a board region. | min (mm) | B | `silk_clearance` against the edge (key `silk_edge_clearance`) | SilkToBoardRegionClearance (value key unknown) | none | KiCad: lossy: shares the silk_clearance value. Altium: lossy: value key unknown. Eagle: none | pcbnew/drc/drc_test_provider_edge_clearance.cpp:176; golden-rules/odrive-v3.txt:10 |
| `text_height` | Height of text, fields, text boxes and dimensions. | min / max (mm) | U (text) | `text_height`; floor `rules.min_text_height` (silk layers) | none | none | KiCad: exact. Altium: none. Eagle: none | common/drc_rules.keywords:60; pcbnew/drc/drc_engine.cpp:209-214 |
| `text_thickness` | Stroke thickness of text. | min / max (mm) | U (text) | `text_thickness`; floor `rules.min_text_thickness` (silk layers) | none | none | KiCad: exact. Altium: none. Eagle: none | common/drc_rules.keywords:61; pcbnew/drc/drc_engine.cpp:216-221 |
| `text_mirroring` | Text mirrored the wrong way for its board side. | variant: front_mirrored, back_unmirrored | U (text) | DRC codes `mirrored_text_on_front_layer`, `nonmirrored_text_on_back_layer` | none | none | KiCad: lossy: board-wide severity only, not scopable. Altium: none. Eagle: none | pcbnew/drc/drc_item.h:127-128 |
| `text_font` | Board text must use a vector (fabrication-safe) font. | none | U (text) | none | none | `checkFont` (flag) | KiCad: none. Altium: none. Eagle: lossy: flag only; meaning unverified (GK) | golden-rules/moteus.txt:91 |

### 3.6 Copper pour, plane and thermal relief

| Canonical id | Meaning | Parameters (units) | Arity | KiCad | Altium | Eagle | Mapping quality | Provenance |
|---|---|---|---|---|---|---|---|---|
| `zone_connection` | How a pad or via joins a pour: solid, thermal relief, none, or thermal for plated through-holes only. | value: solid, thermal, none, thermal_pth_only | B (pad or via, zone) | `zone_connection`; pad, footprint and zone connection settings | PolygonConnect `CONNECTSTYLE`; PlaneConnect `PLANECONNECTSTYLE`; per-object prefixes `THPAD.`, `SMDPAD.`, `VIA.`, `PAD.` | polygon `thermals`; pad `thermals="no"`; `slThermalsForVias` | KiCad: exact (`thermal_pth_only` has no DRU token but exists as an item setting). Altium: exact (prefixes become type() scopes). Eagle: lossy: per polygon and one board-wide via flag only | common/drc_rules.keywords:78; golden-rules/odrive-v3.txt:32, :34; golden-rules/moteus.txt:91 |
| `thermal_relief_gap` | Gap between a pad or via and the pour around a thermal relief. | min (mm) | B (pad or via, zone) | `thermal_relief_gap`; zone thermal gap; pad override | PolygonConnect `AIRGAPWIDTH`; PlaneConnect `RELIEFAIRGAP` | `slThermalIsolate` | KiCad: exact. Altium: exact (PlaneConnect `RELIEFEXPANSION` meaning unverified (GK)). Eagle: exact (meaning from OMX) | common/drc_rules.keywords:63; golden-rules/odrive-v3.txt:32, :34; golden-rules/moteus.txt:91 |
| `thermal_spoke_width` | Width of thermal-relief spokes. | opt, with min / max clamp (mm) | B (pad or via, zone) | `thermal_spoke_width`; zone spoke width; pad override | `RELIEFCONDUCTORWIDTH` (PolygonConnect, PlaneConnect) | none | KiCad: exact. Altium: exact. Eagle: none | common/drc_rules.keywords:64; golden-rules/odrive-v3.txt:32, :34 |
| `thermal_spoke_count` | Generation: number and angle of spokes in a thermal relief. | count (count); angle (deg) | B (pad or via, zone) | pad property `Thermal_Relief_Spoke_Angle` only | `RELIEFENTRIES` (count); PolygonConnect `POLYGONRELIEFANGLE` (deg) | none | KiCad: lossy: angle per pad only, no count. Altium: exact. Eagle: none | golden-rules/odrive-v3.txt:32, :34 |
| `starved_thermal` | A thermal relief ends up with too few connected spokes. | min_resolved (count, 0-4); or min_copper (%) | B (pad, zone) | `min_resolved_spokes`; floor `rules.min_resolved_spokes`; DRC code `starved_thermal` | DRC option `REPORTSTARVEDTHERMALS` with `MINSTARVEDCOPPERPERCENT` (%) | none | KiCad: exact (count form). Altium: lossy: percent threshold; DRC option, not a rule. Eagle: none | common/drc_rules.keywords:40; odrive-DRC D6 |
| `zone_min_width` | Minimum width of filled copper in a pour. | min (mm) | U (zone) | zone minimum width property; `defaults.zones.min_thickness` | none found | polygon `width`; `checkPolygonWidth` (flag) | KiCad: exact (fill parameter, not a DRC check). Altium: none. Eagle: lossy: the pen width also inflates the outline | pcbnew/board_design_settings.cpp:1010-1134; common/io/eagle/eagle_parser.h:968 |
| `zone_fill_style` | Generation settings of a pour: fill mode, hatch, smoothing, island removal, priority, fillets, thieving. | fill_mode (solid, hatch); hatch_width / hatch_gap (mm); hatch_angle (deg); island_removal (always, never, below_area); min_island_area (mm2); priority (count); external_fillets (bool) | U (zone) | zone properties; `defaults.zones.*`; `zones_allow_external_fillets` | none (polygon settings are object properties (GK)) | polygon `pour`, `spacing`, `rank`, `orphans` | KiCad: exact. Altium: none. Eagle: lossy: no smoothing or fillet settings | pcbnew/board_design_settings.cpp:1010-1145; common/io/eagle/eagle_parser.h:961-1004 |
| `via_stitch_overlap` | Same-net via-stitching areas overlap. | none | U | DRC code `via_stitch_overlap` | none | none | KiCad: lossy: board-wide severity only, not scopable. Altium: none. Eagle: none | pcbnew/drc/drc_item.h:139 |

### 3.7 Routing style and topology

| Canonical id | Meaning | Parameters (units) | Arity | KiCad | Altium | Eagle | Mapping quality | Provenance |
|---|---|---|---|---|---|---|---|---|
| `routing_layers` | Copper layers on which the scoped nets may be routed. | allowed_layers (list) | U (net) | `disallow track` on the other layers | RoutingLayers (one bool per layer) | none | KiCad: lossy: written as a prohibition on the complement; router intent lost. Altium: exact (mid-layer names need the layer stack). Eagle: none | golden-rules/odrive-v3.txt:29; common/drc_rules.keywords:16 |
| `routing_topology` | Required net topology. | value: shortest, horizontal, vertical, daisychain_simple, daisychain_middriven, daisychain_balanced, starburst | U (net) | none | RoutingTopology (`TOPOLOGY`) | none | KiCad: none. Altium: exact. Eagle: none | golden-rules/odrive-v3.txt:24 |
| `routing_priority` | Order in which the autorouter treats nets. | value (count) | U (net) | none | RoutingPriority (`ROUTINGPRIORITY`) | none | KiCad: none. Altium: exact. Eagle: none | golden-rules/odrive-v3.txt:30 |
| `routing_corners` | Corner style and setback of routed tracks. | style (90, 45, round); setback min / max (mm) | U (track) | none (export approximation: `track_angle`) | RoutingCorners (`CORNERSTYLE`, `MINSETBACK`, `MAXSETBACK`) | `checkAngle` (flag) | KiCad: none. Altium: exact. Eagle: lossy: flag for off-45 deg wires only; meaning unverified (GK) | golden-rules/odrive-v3.txt:28; golden-rules/moteus.txt:91 |
| `fanout_control` | Fan-out pattern for component pads. | style (auto, rows, staggered, bga, underpads); direction (enum); bga_dir, bga_via_mode (string); via_grid (mm) | U (component) | none | FanoutControl (`FANOUTSTYLE`, `FANOUTDIRECTION`, `BGADIR`, `BGAVIAMODE`, `VIAGRID`) | none | KiCad: none. Altium: exact. Eagle: none | golden-rules/odrive-v3.txt:16-19 |
| `smd_entry` | Directions from which a track may enter an SMD pad. | side, corner, any_angle (bool) | U (pad) | none | SMDEntry (`SIDE`, `CORNER`, `ANYANGLE`) | none | KiCad: none. Altium: exact. Eagle: none | W6:4; AD10 |
| `smd_to_corner` | Minimum distance from an SMD pad edge to the first track corner. | min (mm) | U (pad) | none | SMDToCorner (`DISTANCE`, stored without unit) | none | KiCad: none. Altium: lossy: unit of `DISTANCE` inferred. Eagle: none | W6:24; AD10 |
| `smd_to_plane` | Maximum routed length from an SMD pad to its plane connection. | max (mm) | U (pad) | none | SMDToPlane (key unknown) | none | KiCad: none. Altium: lossy: key unknown. Eagle: none | AD10 |
| `smd_neckdown` | Maximum ratio of track width to SMD pad width. | max (%) | U (pad) | none | SMDNeckDown (key unknown) | none | KiCad: none. Altium: lossy: key unknown. Eagle: none | AD10 |
| `routing_neckdown` | Altium Routing Neck-Down; meaning unknown. | unknown | unknown | none | RoutingNeckDown (no record or documentation found) | none | KiCad: none. Altium: lossy: meaning unknown. Eagle: none | unverified (L1: confirm via Altium scripting, Q-A12) |
| `track_not_centered_on_via` | A track ends near, but not on, a via centre. | none | B (track, via) | DRC code `track_not_centered_on_via` | none | none | KiCad: lossy: board-wide severity only, not scopable. Altium: none. Eagle: none | pcbnew/drc/drc_item.h:135 |
| `router_presets` | Preferred widths, via sizes and pair geometries offered to the router; not checked. | lists of width (mm); via diameter / drill (mm); pair width / gap / via gap (mm) | - | `track_widths`, `via_dimensions`, `diff_pair_dimensions`, `via_stack_presets` | none (preferred values live inside rules) | none | KiCad: exact. Altium: none. Eagle: none | pcbnew/board_design_settings.cpp:432-653 |
| `teardrop_style` | Teardrop generator settings; not checked. | target flags; max length / height (mm); ratios | - | `teardrop_options`, `teardrop_parameters`; item teardrop properties | none found | none | KiCad: exact. Altium: none. Eagle: none | pcbnew/board_design_settings.cpp:656-798 |

### 3.8 High speed: length, differential pair, skew

| Canonical id | Meaning | Parameters (units) | Arity | KiCad | Altium | Eagle | Mapping quality | Provenance |
|---|---|---|---|---|---|---|---|---|
| `length` | Routed length (or delay) of each scoped net. | min / opt / max (mm or ps); opt = tuning target; include_via_height (bool) | U (net) | `length`; `rules.use_height_for_length_calcs` | Length (keys unknown) | none | KiCad: exact. Altium: lossy: keys unknown. Eagle: none | common/drc_rules.keywords:30; AD4 |
| `net_chain_length` | Trunk length of a net chain between its terminal pads, through series parts. | min / max (mm or ps) | U (chain) | `net_chain_length` (master only) | none | none | KiCad: exact. Altium: none. Eagle: none | common/drc_rules.keywords:31 |
| `stub_length` | Length of stubs off a chain (KiCad: stub member nets; Altium: daisy-chain stubs). | max (mm or ps) | U | `stub_length` (master only) | StubLength (key unknown) | none | KiCad: exact. Altium: lossy: key unknown; stub definitions differ. Eagle: none | common/drc_rules.keywords:33; AD8 |
| `return_path` | Reference plane continues under the routed net. | reference_layer; reference_net (pattern); min_margin (mm); exclude_voids (bool) | U (net) | `return_path` (layer, net; coverage only; master only) | ReturnPath (keys unknown) | none | KiCad: lossy: no margin and no void option. Altium: lossy: keys unknown; reference layer comes from the impedance profile. Eagle: none | common/drc_rules.keywords:32; AD6 |
| `skew` | Length or delay difference among the scoped nets, or within each differential pair. | max, or min / opt / max (mm or ps); within_diff_pairs (bool); compare (set: in_pair, pair_vs_pair, xsignals, others) | U (net group) | `skew` with `within_diff_pairs` | MatchedLengths (`TOLERANCE`, `CHECK*` flags) | `dpMaxLengthDifference` | KiCad: lossy: no compare selection. Altium: exact. Eagle: lossy: within pairs only; meaning unverified (GK) | common/drc_rules.keywords:53; W5:7-8; golden-rules/moteus.txt:91 |
| `diff_pair_gap` | Gap between the two tracks of a coupled differential pair. | min / opt / max (mm); via_gap (mm, router only) | B (pair) | `diff_pair_gap`; netclass `diff_pair_gap` (opt), `diff_pair_via_gap` | DiffPairsRouting (`MINLIMIT`, `MAXLIMIT`, `MOSTFREQGAP`) | `dpGapFactor` | KiCad: exact. Altium: exact (gap reading from third-party). Eagle: lossy: base of the factor unknown | common/drc_rules.keywords:14; golden-rules/odrive-v3.txt:15; golden-rules/moteus.txt:91 |
| `diff_pair_uncoupled` | Length a differential-pair track runs uncoupled from its partner. | max (mm) | B (pair) | `diff_pair_uncoupled` | DiffPairsRouting `MAXUNCOUPLEDLENGTH` | none | KiCad: exact. Altium: exact. Eagle: none | common/drc_rules.keywords:15; golden/odrive-v3/src/v3/PCB.PcbDoc rec DiffPairsRouting |
| `parallel_segment` | Summed length of segments running in parallel within a gap (crosstalk guard). | gap (mm); max_parallel_length (mm); layer_mode (same, adjacent) | B | none | ParallelSegment (keys unknown) | none | KiCad: none. Altium: lossy: keys unknown. Eagle: none | AD8 |
| `via_count` | Number of vias on each scoped net. | min / max (count) | U (net) | `via_count` | MaximumViaCount (max; key unknown) | none | KiCad: exact. Altium: lossy: max only; key unknown. Eagle: none | common/drc_rules.keywords:72; AD4 |
| `backdrill_stub` | Maximum via stub length; longer stubs make the via a back-drill candidate. | max_stub (mm); oversize (mm) | U (via) | none (back-drill properties exist, no check) | BackDrilling (keys unknown) | none | KiCad: none. Altium: lossy: keys unknown. Eagle: none | AD7 |
| `impedance_target` | Controlled-impedance target of routed geometry, single-ended or differential. | min / opt / max (ohm); profile (id) | U (net) | netclass `tuning_profile` and project tuning profiles (derived width and gap per layer) | Width `MINIMP`, `MAXIMP`, `FAVIMP`, `IMPEDANCEPROFILE*`; DiffPairsRouting `IMPEDANCEPROFILE*` | none | KiCad: lossy: KiCad keeps the derived geometry, not the ohm target. Altium: exact (profile needs the layer stack). Eagle: none | pcbnew/drc/drc_engine.cpp:473-598; golden/odrive-v3/src/v3/PCB.PcbDoc rec Width; W5:2, :32 |
| `tuning_profile_consistency` | Netclass and net-chain tuning profiles are defined and consistent. | variant: missing_profile, chain_mixed_profiles, implicit_geometry | U | DRC codes `missing_tuning_profile`, `net_chain_tuning_profiles`; key `tuning_profile_track_geometries` | none | none | KiCad: exact. Altium: none. Eagle: none | pcbnew/drc/drc_item.h:121, :130-131 |
| `meander_style` | Length-tuning (meander) generator settings; not checked. | min / max amplitude (mm); spacing (mm); corner_radius (%); single_sided (bool) | - | `tuning_pattern_settings`; tuning-pattern item properties | none found | none | KiCad: exact. Altium: none. Eagle: none | pcbnew/board_design_settings.cpp:800-866 |

### 3.9 Placement, courtyard, height and keepout

| Canonical id | Meaning | Parameters (units) | Arity | KiCad | Altium | Eagle | Mapping quality | Provenance |
|---|---|---|---|---|---|---|---|---|
| `courtyard_clearance` | Minimum clearance between component outlines (KiCad: 2D courtyards per side; Altium: component bodies, horizontal and vertical). | min (mm); vertical_min (mm); check_mode (enum) | B (footprint pair) | `courtyard_clearance` (key `courtyards_overlap`) | ComponentClearance (`GAP`, `VERTICALGAP`, `COLLISIONCHECKMODE`, `SHOWDISTANCES`) | none | KiCad: exact (2D form). Altium: lossy: Altium measures 3D bodies; `COLLISIONCHECKMODE=3` unresolved. Eagle: none | common/drc_rules.keywords:12; golden-rules/odrive-v3.txt:23, :38-39 |
| `courtyard_required` | Every footprint has a valid courtyard. | variant: missing, malformed | U (footprint) | DRC codes `missing_courtyard`, `malformed_courtyard`; footprint courtyard exemption | none | none | KiCad: lossy: board-wide severity only, not scopable. Altium: none. Eagle: none | pcbnew/drc/drc_item.h:80-81 |
| `hole_in_courtyard` | A plated or unplated hole sits inside another footprint's courtyard. | variant: pth, npth | B (hole, footprint) | DRC codes `pth_inside_courtyard`, `npth_inside_courtyard` | none | none | KiCad: lossy: board-wide severity only, not scopable. Altium: none. Eagle: none | pcbnew/drc/drc_item.h:83-84 |
| `component_height` | Height of components, usually scoped to an area. | min / opt / max (mm) | U (component) | none as a check (expressible as an `assertion` on footprint `Height`) | Height (`MINHEIGHT`, `PREFHEIGHT`, `MAXHEIGHT`) | none | KiCad: lossy: assertion form, no opt. Altium: exact. Eagle: none | golden-rules/odrive-v3.txt:22; AD11 |
| `room_confinement` | Scoped components must stay inside (or outside) a region. Every matching rule applies (all-match, section 5, P12). | region (polygon, mm); mode (confine_in, confine_out); lock (bool) | U (component) | rule area plus `assertion` (`enclosedByArea`) or `disallow footprint` | RoomDefinition (vertices, `CONFINEMENTSTYLE`, `LOCKCOMPONENTS`) | none | KiCad: lossy: area plus assertion. Altium: exact. Eagle: none | W8 rec RoomDefinition; AD11 |
| `component_orientation` | Allowed rotations of scoped components. | allowed (set of deg) | U (component) | `assertion` on `Orientation` (graphical editor output) | ComponentOrientations (keys unknown) | none | KiCad: lossy: assertion form. Altium: lossy: keys unknown. Eagle: none | pcbnew/drc/rule_editor/drc_re_allowed_orientation_constraint_data.h:95; AD4 |
| `permitted_layers` | Board sides on which scoped components may be placed. | allowed (set: top, bottom) | U (component) | `assertion` (graphical editor output) | PermittedLayers (keys unknown) | none | KiCad: lossy: assertion form. Altium: lossy: keys unknown. Eagle: none | pcbnew/drc/rule_editor/drc_re_permitted_layers_constraint_data.h:76; AD4 |
| `nets_to_ignore` | Nets that placement tools ignore; no DRC check. | none (scope only) | U (net) | none | NetsToIgnore (keys unknown) | none | KiCad: none. Altium: lossy: meaning unverified (GK). Eagle: none | AD4 |
| `disallow` | Object kinds forbidden in a region, on a layer or in a scope (keepout). | kinds (set: track, via, through_via, blind_via, buried_via, micro_via, pad, zone_fill, text, graphic, hole, footprint) | U | `disallow`; keepout rule areas | none as a rule (keepouts are primitives with the `KEEPOUT` flag; verify in L1) | tRestrict, bRestrict, vRestrict shapes; polygon `pour="cutout"`; `checkRestrict` | KiCad: exact. Altium: none. Eagle: lossy: restrict layers forbid copper or vias only | common/drc_rules.keywords:16; pcbnew/drc/drc_engine.cpp:600-678; pcbnew/pcb_io/eagle/pcb_io_eagle.cpp:121-147 |

### 3.10 Testpoint and assembly

| Canonical id | Meaning | Parameters (units) | Arity | KiCad | Altium | Eagle | Mapping quality | Provenance |
|---|---|---|---|---|---|---|---|---|
| `testpoint_style` | Allowed testpoint pad and hole sizes, grid and sides. | variant: fabrication, assembly; size min / opt / max (mm); hole min / opt / max (mm); grid, grid_tolerance (mm); sides (set); under_component (bool) | U (pad) | none (pads carry a testpoint fabrication property, no rule) | FabricationTestpoint, AssemblyTestpoint | none | KiCad: none. Altium: lossy: `SIDE=3` unresolved. Eagle: none | golden-rules/odrive-v3.txt:4-5 |
| `testpoint_usage` | Whether each scoped net needs a testpoint, and whether several are allowed. | variant: fabrication, assembly; value (required, invalid, ignore); allow_multiple (bool) | U (net) | none | FabricationTestPointUsage, AssemblyTestPointUsage (`VALID`, `ALLOWMULTIPLE`) | none | KiCad: none. Altium: lossy: `VALID` codes read from enum order only. Eagle: none | golden-rules/odrive-v3.txt:6-7 |
| `part_value_present` | Every placed part has a value (assembly data). | none | U (component) | none | none | `checkValues` (flag) | KiCad: none. Altium: none. Eagle: lossy: flag only; meaning unverified (GK) | golden-rules/moteus.txt:91 |

### 3.11 Manufacturing (fab capability)

Most fab capabilities are floors on kinds defined elsewhere (track width, clearance, hole size, annular width, mask
sliver, hole-to-hole, edge clearance, aspect ratio); section 7.1 shows how they attach. The kinds below have no
other home. `silk_line_width` and `board_size` are project-defined for W7: no KiCad, Altium or Eagle rule defines
them in the sources read.

| Canonical id | Meaning | Parameters (units) | Arity | KiCad | Altium | Eagle | Mapping quality | Provenance |
|---|---|---|---|---|---|---|---|---|
| `outline_integrity` | The board outline is valid and objects respect it. | variant: invalid_outline, copper_outside_board, text_on_edge_cuts | U | DRC codes `invalid_outline`, `copper_outside_board`, `text_on_edge_cuts` | none found | none | KiCad: lossy: board-wide severity only, not scopable. Altium: none. Eagle: none | pcbnew/drc/drc_item.h:51, :86, :142 |
| `copper_sliver` | Narrow acute copper slivers. | width_tolerance (mm); angle_tolerance (deg); min_length (mm) | U | DRC code `copper_sliver` (thresholds in advanced config) | none | none | KiCad: lossy: thresholds not stored in the project. Altium: none. Eagle: none | pcbnew/drc/drc_test_provider_sliver_checker.cpp:88, :91, :193 |
| `silk_line_width` | Minimum silkscreen line width a fab can print. Project-defined for W7; no tool rule found. | min (mm) | U (graphic) | none | none | none | KiCad: none. Altium: none. Eagle: none | pcb-second-brain/docs/02-knowledge-base.md §2.5 |
| `board_size` | Board outline must fit the fab's size limits. Project-defined for W7; no tool rule found. | max_x / max_y, min_x / min_y (mm) | U (board) | none | none | none | KiCad: none. Altium: none. Eagle: none | pcb-second-brain/docs/02-knowledge-base.md §2.5 |
| `stackup` | The stackup (layer count, copper and dielectric thickness, materials) must be one the fab offers. | layers (count); per-layer thickness (mm) and material; board thickness (mm) | U (board) | board `(setup (stackup))`; KiCad 5 `.pro` `CopperLayerCount`, `BoardThickness` | LayerStack (legacy kind, keys unknown); `Board6/Data` stack (not decoded) | `layerSetup`, `mtCopper`, `mtIsolate` | KiCad: exact. Altium: lossy: not decoded yet. Eagle: lossy: syntax and indexing unverified (GK) | include/board_design_settings.h:782; W1:54; golden-rules/moteus.txt:91 |

### 3.12 Signal integrity

All from Altium; KiCad and Eagle have no signal-integrity rules. These are checked by simulation (T4), not by DRC.

| Canonical id | Meaning | Parameters (units) | Arity | KiCad | Altium | Eagle | Mapping quality | Provenance |
|---|---|---|---|---|---|---|---|---|
| `si_stimulus` | Stimulus used by signal-integrity analysis. | kind (constant_level, single_pulse, periodic_pulse); levels (V); start / stop / period (s) | U (net) | none | SignalStimulus | none | KiCad: none. Altium: lossy: keys unknown (AD9 summary only). Eagle: none | AD9 |
| `si_overshoot` | Maximum overshoot on an edge. | variant: rising, falling; max (V) | U (net) | none | OvershootRising, OvershootFalling | none | KiCad: none. Altium: lossy: keys unknown (AD9 summary only). Eagle: none | AD9 |
| `si_undershoot` | Maximum undershoot on an edge. | variant: rising, falling; max (V) | U (net) | none | UndershootRising, UndershootFalling | none | KiCad: none. Altium: lossy: keys unknown (AD9 summary only). Eagle: none | AD9 |
| `si_impedance` | Net impedance from signal-integrity analysis (unlike impedance_target, a simulated result). | min / max (ohm; unit unverified (GK)) | U (net) | none | MaxMinImpedance | none | KiCad: none. Altium: lossy: keys unknown (AD9 summary only). Eagle: none | AD9 |
| `si_signal_top_value` | Minimum settled high level. | min (V) | U (net) | none | SignalTopValue | none | KiCad: none. Altium: lossy: keys unknown (AD9 summary only). Eagle: none | AD9 |
| `si_signal_base_value` | Maximum settled low level. | max (V) | U (net) | none | SignalBaseValue | none | KiCad: none. Altium: lossy: keys unknown (AD9 summary only). Eagle: none | AD9 |
| `si_flight_time` | Maximum interconnect delay on an edge. | variant: rising, falling; max (s) | U (net) | none | FlightTimeRising, FlightTimeFalling | none | KiCad: none. Altium: lossy: keys unknown (AD9 summary only). Eagle: none | AD9 |
| `si_slope` | Maximum time from threshold to valid level on an edge. | variant: rising, falling; max (s) | U (net) | none | SlopeRising, SlopeFalling | none | KiCad: none. Altium: lossy: keys unknown (AD9 summary only). Eagle: none | AD9 |
| `supply_net_voltage` | Nominal voltage of a supply net, for SI and power checks. | value (V) | U (net) | none (no supply-voltage attribute in DRU) | SupplyNets (`VOLTAGE`) | none | KiCad: none. Altium: exact. Eagle: none | W9 rec SupplyNets; AD9 |

### 3.13 Schematic and ERC

KiCad's 63 ERC codes group into these kinds; the per-code table is section 4.1.13. Altium's ERC levels are stored per
violation type number, and no source gives the type names (Q-A6), so every Altium entry here is lossy until L1.

| Canonical id | Meaning | Parameters (units) | Arity | KiCad | Altium | Eagle | Mapping quality | Provenance |
|---|---|---|---|---|---|---|---|---|
| `erc_pin_type_conflict` | Two connected pins (or ports) have conflicting electrical types. | matrix (type x type: ok, warning, error) | B (pin pair) | `erc.pin_map` (12 x 12); key `pin_to_pin` (warning and error codes) | `[ERC Connection Matrix]` (17 x 17) | none | KiCad: exact. Altium: lossy: 17 types incl. ports and sheet entries; labels unverified (GK). Eagle: none | eeschema/erc/erc_settings.cpp:43-58; golden/odrive-v3/src/v3/two_ax_PCB.PrjPcb:4620-4637 |
| `erc_undriven_input` | An input or power-input pin has no driver on its net. | variant: input, power_input | U (net) | `pin_not_driven`, `power_pin_not_driven` | `[Electrical Rules Check]` levels | none | KiCad: exact. Altium: lossy: Altium `Type<k>` numbers not mapped (Q-A6). Eagle: none | eeschema/erc/erc_settings.h:46, :48; eeschema/erc/erc.cpp:307-332 |
| `erc_unconnected_pin` | A pin is unconnected and carries no no-connect flag. | none | U (pin) | `pin_not_connected` | `[Electrical Rules Check]` levels | none | KiCad: exact. Altium: lossy: Altium `Type<k>` numbers not mapped (Q-A6). Eagle: none | eeschema/erc/erc_settings.h:45 |
| `erc_no_connect` | No-connect flag misuse. | variant: connected, dangling | U | `no_connect_connected`, `no_connect_dangling` | `[Electrical Rules Check]` levels | none | KiCad: exact. Altium: lossy: Altium `Type<k>` numbers not mapped (Q-A6). Eagle: none | eeschema/erc/erc_settings.h:51-52 |
| `erc_dangling` | Dangling schematic items. | variant: label, wire, wire_endpoint | U | `label_dangling`, `wire_dangling`, `unconnected_wire_endpoint` | `[Electrical Rules Check]` levels | none | KiCad: exact. Altium: lossy: Altium `Type<k>` numbers not mapped (Q-A6). Eagle: none | eeschema/erc/erc_settings.h:53, :81, :94 |
| `erc_label` | Label naming and label connectivity problems. | variant: similar_labels, similar_power, similar_label_and_power, single_global_label, same_local_global_label, same_local_global_power, isolated_pin_label, label_multiple_wires, empty_label_name, multiple_net_names, hier_label_mismatch | U | 11 ERC keys (section 4.1.13) | `[Electrical Rules Check]` levels | none | KiCad: exact. Altium: lossy: Altium `Type<k>` numbers not mapped (Q-A6). Eagle: none | eeschema/erc/erc_settings.h:50-59, :70, :77, :93, :101 |
| `erc_bus` | Bus definition and bus connection problems. | variant: definition_conflict, net_not_bus_member, bus_to_bus, bus_to_net, entry_needed | U | `bus_definition_conflict`, `net_not_bus_member`, `bus_to_bus_conflict`, `bus_to_net_conflict`, `bus_entry_needed` | `[Electrical Rules Check]` levels | none | KiCad: exact. Altium: lossy: Altium `Type<k>` numbers not mapped (Q-A6). Eagle: none | eeschema/erc/erc_settings.h:69-74, :91 |
| `erc_annotation` | Reference designators and sheet instances are complete and unique. | variant: unannotated, duplicate_reference, extra_units, unit_value_mismatch, duplicate_sheet_name, duplicate_sheet_suffix, ref_suffix_mismatch, missing_name | U | 7 ERC keys (section 4.1.13) | `[Electrical Rules Check]` levels | `checkNames` (board-side flag) | KiCad: exact. Altium: lossy: Altium `Type<k>` numbers not mapped (Q-A6). Eagle: lossy: board-side on/off flag only; meaning unverified (GK) | eeschema/erc/erc_settings.h:43, :87-90, :106-107; golden-rules/moteus.txt:91 |
| `erc_multi_unit` | Units of multi-unit symbols are consistent and placed. | variant: different_unit_footprint, different_unit_net, missing_unit, missing_power_pin, missing_input_pin, missing_bidi_pin, duplicate_pins | U | 7 ERC keys (section 4.1.13) | `[Electrical Rules Check]` levels | none | KiCad: exact. Altium: lossy: Altium `Type<k>` numbers not mapped (Q-A6). Eagle: none | eeschema/erc/erc_settings.h:60-67, :111 |
| `erc_symbol_library` | Symbol, footprint-link and pin-map integrity. | variant: lib_symbol_issues, lib_symbol_mismatch, footprint_link_issues, footprint_filter, variant_symbol_invalid, variant_symbol_incompatible, stacked_pin_name, field_name_whitespace, pin_map_bad_pad, pin_map_unmapped_pin, pin_map_duplicate_pad, pin_map_stale_pin | U | 12 ERC keys (section 4.1.13) | `[Electrical Rules Check]` levels | none | KiCad: exact. Altium: lossy: Altium `Type<k>` numbers not mapped (Q-A6). Eagle: none | eeschema/erc/erc_settings.h:82-86, :95-103 |
| `erc_power` | Power-net specific problems. | variant: ground_pin_not_ground, wired_implicit_power, empty_power_value | U | `ground_pin_not_ground`, `wired_implicit_power`, `empty_power_value` | `[Electrical Rules Check]` levels | none | KiCad: exact. Altium: lossy: Altium `Type<k>` numbers not mapped (Q-A6). Eagle: none | eeschema/erc/erc_settings.h:76, :104-105 |
| `erc_netclass` | A netclass is referenced but not defined. | none | U | `undefined_netclass` | `[Electrical Rules Check]` levels | none | KiCad: exact. Altium: lossy: Altium `Type<k>` numbers not mapped (Q-A6). Eagle: none | eeschema/erc/erc_settings.h:79 |
| `erc_drawing` | Drawing hygiene: off-grid endpoints, four-way junctions. | variant: endpoint_off_grid, four_way_junction | U | `endpoint_off_grid`, `four_way_junction` | `[Electrical Rules Check]` levels | none | KiCad: exact. Altium: lossy: Altium `Type<k>` numbers not mapped (Q-A6). Eagle: none | eeschema/erc/erc_settings.h:44, :92 |
| `erc_simulation_model` | A symbol's simulation model has a problem. | none | U (symbol) | `simulation_model_issue` | `[Electrical Rules Check]` levels | none | KiCad: exact. Altium: lossy: Altium `Type<k>` numbers not mapped (Q-A6). Eagle: none | eeschema/erc/erc_settings.h:80 |

### 3.14 General: custom assertions, data integrity, engine settings

| Canonical id | Meaning | Parameters (units) | Arity | KiCad | Altium | Eagle | Mapping quality | Provenance |
|---|---|---|---|---|---|---|---|---|
| `custom_assertion` | A user expression that every matching object must satisfy. | expr (canonical expression) | U | `assertion`; DRC `assertion_failure`, `generic_warning`, `generic_error`; ERC `generic-warning`, `generic-error` | none found | none | KiCad: exact. Altium: none. Eagle: none | common/drc_rules.keywords:2; pcbnew/drc/drc_engine.cpp:2278-2347 |
| `text_variable_unresolved` | A text variable reference cannot be resolved. | none | U (text) | DRC and ERC `unresolved_variable` | none | none | KiCad: lossy: board-wide severity only, not scopable. Altium: none. Eagle: none | pcbnew/drc/drc_item.h:102; eeschema/erc/erc_settings.h:78 |
| `item_on_disabled_layer` | An object sits on a layer that is not enabled in the stackup. | none | U | DRC code `item_on_disabled_layer` | none | none | KiCad: lossy: board-wide severity only, not scopable. Altium: none. Eagle: none | pcbnew/drc/drc_item.h:85 |
| `footprint_integrity` | Footprint and padstack sanity. | variant: type_mismatch, th_pad_without_hole, padstack, padstack_invalid, net_tie_definition, scaled_with_pads, drilled_smd_pad, invalid_multilayer_pad | U (footprint) | 6 DRC codes of the `footprint checks` provider | DRC options `REPORTDRILLEDSMTPADS`, `REPORTINVALIDMULTILAYERPADS` | none | KiCad: lossy: board-wide severity only, not scopable. Altium: lossy: DRC options; meanings unverified (GK). Eagle: none | pcbnew/drc/drc_test_provider_footprint_checks.cpp:28-39; odrive-DRC D4 |
| `library_parity` | A footprint is missing from its library or differs from the library copy. | variant: missing, mismatch | U (footprint) | DRC codes `lib_footprint_issues`, `lib_footprint_mismatch` | none | none | KiCad: lossy: board-wide severity only, not scopable. Altium: none. Eagle: none | pcbnew/drc/drc_item.h:96-97 |
| `board_off_grid` | Board objects lie off the placement grid. | grid (mm) | U | none | none | `checkGrid` (flag) | KiCad: none. Altium: none. Eagle: lossy: flag only; meaning unverified (GK) | golden-rules/moteus.txt:103 |
| `drc_engine_settings` | Not a rule: DRC accuracy and run settings carried with a rule set. | arc_tolerance (mm); epsilon (mm); hole_wall_thickness (mm); max_violations (count); checks_enabled (set) | - | `rules.max_error`; advanced config epsilon and hole wall; `rule_severities` | `MAXVIOLATIONCOUNT`, `RULESETTOCHECK`, `ONLINERULESETTOCHECK`, report options | `maxErrors`; `check*` switches | KiCad: exact (as rule-set metadata). Altium: exact (as rule-set metadata). Eagle: exact (as rule-set metadata) | pcbnew/board_design_settings.cpp:221, :2054-2063; odrive-DRC D1-D3, D8; golden-rules/moteus.txt:91 |

## 4. Bidirectional mapping tables

Each table lists every item of its kind found in the notes; nothing is left out. The "Canonical target" column is
the import direction; section 3 gives the export quality per kind, and section 4.4 the export recipes.

### 4.1 KiCad -> canonical

#### 4.1.1 DRU constraint types (all 42 `DRC_CONSTRAINT_T` values and the 3 deprecated spellings)

| # | DRU keyword | Enum | Canonical target | Quality | Note | Provenance |
|---|---|---|---|---|---|---|
| 1 | `clearance` | `CLEARANCE_CONSTRAINT` | clearance.min | exact | a rule with min <= 0 mm or `(severity ignore)` also implies short_circuit.allowed=true and no tracks_crossing check for its scope (section 3.1) | common/drc_rules.keywords:8; pcbnew/drc/drc_constraint_types.h:26 |
| 2 | `creepage` | `CREEPAGE_CONSTRAINT` | creepage.min | exact |  | common/drc_rules.keywords:13; pcbnew/drc/drc_constraint_types.h:27 |
| 3 | `hole_clearance` | `HOLE_CLEARANCE_CONSTRAINT` | hole_clearance.min | exact |  | common/drc_rules.keywords:24; pcbnew/drc/drc_constraint_types.h:28 |
| 4 | `hole_to_hole` | `HOLE_TO_HOLE_CONSTRAINT` | hole_to_hole.min | exact |  | common/drc_rules.keywords:26; pcbnew/drc/drc_constraint_types.h:29 |
| 5 | `edge_clearance` | `EDGE_CLEARANCE_CONSTRAINT` | edge_clearance.min | exact |  | common/drc_rules.keywords:17; pcbnew/drc/drc_constraint_types.h:30 |
| 6 | `physical_clearance` | `PHYSICAL_CLEARANCE_CONSTRAINT` | physical_clearance.min | exact |  | common/drc_rules.keywords:48; pcbnew/drc/drc_constraint_types.h:57 |
| 7 | `physical_hole_clearance` | `PHYSICAL_HOLE_CLEARANCE_CONSTRAINT` | physical_hole_clearance.min | exact |  | common/drc_rules.keywords:49; pcbnew/drc/drc_constraint_types.h:58 |
| 8 | `courtyard_clearance` | `COURTYARD_CLEARANCE_CONSTRAINT` | courtyard_clearance.min | exact |  | common/drc_rules.keywords:12; pcbnew/drc/drc_constraint_types.h:32 |
| 9 | `silk_clearance` | `SILK_CLEARANCE_CONSTRAINT` | silk_clearance.min + silk_to_mask_clearance.min + silk_to_edge_clearance.min | exact | split by the target object of B: silk, mask or copper, edge | common/drc_rules.keywords:52; pcbnew/drc/drc_constraint_types.h:33 |
| 10 | `hole_size` | `HOLE_SIZE_CONSTRAINT` | hole_size.min/opt/max | exact |  | common/drc_rules.keywords:25; pcbnew/drc/drc_constraint_types.h:31 |
| 11 | `annular_width` | `ANNULAR_WIDTH_CONSTRAINT` | annular_width.min/max | exact |  | common/drc_rules.keywords:1; pcbnew/drc/drc_constraint_types.h:38 |
| 12 | `via_diameter` | `VIA_DIAMETER_CONSTRAINT` | via_diameter.min/opt/max | exact |  | common/drc_rules.keywords:74; pcbnew/drc/drc_constraint_types.h:47 |
| 13 | `track_width` | `TRACK_WIDTH_CONSTRAINT` | track_width.min/opt/max | exact |  | common/drc_rules.keywords:68; pcbnew/drc/drc_constraint_types.h:36 |
| 14 | `track_segment_length` | `TRACK_SEGMENT_LENGTH_CONSTRAINT` | track_segment_length.min/max | exact |  | common/drc_rules.keywords:69; pcbnew/drc/drc_constraint_types.h:37 |
| 15 | `track_angle` | `TRACK_ANGLE_CONSTRAINT` | track_angle.min/max (deg) | exact |  | common/drc_rules.keywords:67; pcbnew/drc/drc_constraint_types.h:61 |
| 16 | `connection_width` | `CONNECTION_WIDTH_CONSTRAINT` | connection_width.min | exact |  | common/drc_rules.keywords:10; pcbnew/drc/drc_constraint_types.h:60 |
| 17 | `text_height` | `TEXT_HEIGHT_CONSTRAINT` | text_height.min/max | exact |  | common/drc_rules.keywords:60; pcbnew/drc/drc_constraint_types.h:34 |
| 18 | `text_thickness` | `TEXT_THICKNESS_CONSTRAINT` | text_thickness.min/max | exact |  | common/drc_rules.keywords:61; pcbnew/drc/drc_constraint_types.h:35 |
| 19 | `zone_connection` | `ZONE_CONNECTION_CONSTRAINT` | zone_connection.value | exact | solid, thermal_reliefs, none map to solid, thermal, none | common/drc_rules.keywords:78; common/drc_rules.keywords:58 |
| 20 | `thermal_relief_gap` | `THERMAL_RELIEF_GAP_CONSTRAINT` | thermal_relief_gap.min | exact |  | common/drc_rules.keywords:63; pcbnew/drc/drc_constraint_types.h:40 |
| 21 | `thermal_spoke_width` | `THERMAL_SPOKE_WIDTH_CONSTRAINT` | thermal_spoke_width.opt (min/max clamp) | exact |  | common/drc_rules.keywords:64; pcbnew/drc/drc_constraint_types.h:41 |
| 22 | `min_resolved_spokes` | `MIN_RESOLVED_SPOKES_CONSTRAINT` | starved_thermal.min_resolved (count) | exact |  | common/drc_rules.keywords:40; pcbnew/drc/drc_constraint_types.h:42 |
| 23 | `solder_mask_expansion` | `SOLDER_MASK_EXPANSION_CONSTRAINT` | solder_mask_expansion.opt | exact |  | common/drc_rules.keywords:54; pcbnew/drc/drc_constraint_types.h:43 |
| 24 | `solder_paste_abs_margin` | `SOLDER_PASTE_ABS_MARGIN_CONSTRAINT` | solder_paste_margin.abs | exact |  | common/drc_rules.keywords:55; pcbnew/drc/drc_constraint_types.h:44 |
| 25 | `solder_paste_rel_margin` | `SOLDER_PASTE_REL_MARGIN_CONSTRAINT` | solder_paste_margin.ratio | lossy: written value encoding unresolved (Q-K1) |  | common/drc_rules.keywords:56; pcbnew/drc/drc_constraint_types.h:45 |
| 26 | `bridged_mask` | `BRIDGED_MASK_CONSTRAINT` | with `(severity ignore)`: solder_mask_bridge.allowed=true, scope.expr = scope.expr_b = the condition (both participants must match), options.aperture_alone=true (a bare mask aperture matching the condition is exempt on its own). With any other severity, or none: the rule exempts nothing in KiCad, so it is imported with allowed=false and enabled=false, keeping its severity and text in native | lossy: KiCad's constraint is unary and is evaluated once per participant with B = null; only `ignore` is tested; bridges are reported against the provider's internal rule, not the bridged_mask rule | the canonical kind is binary (section 3.5) | common/drc_rules.keywords:5; pcbnew/drc/drc_constraint_types.h:63; pcbnew/drc/drc_rule.cpp:111; pcbnew/drc/drc_test_provider_solder_mask.cpp:379-390, :594-610, :840 |
| 27 | `solder_mask_sliver` | `SOLDER_MASK_SLIVER_CONSTRAINT` | solder_mask_sliver.min | lossy: KiCad parses it but no consumer enforces it (Q-K2) |  | common/drc_rules.keywords:57; pcbnew/drc/drc_constraint_types.h:64 |
| 28 | `disallow` | `DISALLOW_CONSTRAINT` | disallow.kinds (per-token table in section 4.1.2); editor patterns recognised on import: vias under SMD -> vias_under_smd, tracks off permitted layers -> routing_layers, footprints kept out of an area -> room_confinement | exact |  | common/drc_rules.keywords:16; pcbnew/drc/drc_constraint_types.h:46 |
| 29 | `length` | `LENGTH_CONSTRAINT` | length.min/opt/max (mm or ps) | exact |  | common/drc_rules.keywords:30; pcbnew/drc/drc_constraint_types.h:48 |
| 30 | `net_chain_length` | `NET_CHAIN_LENGTH_CONSTRAINT` | net_chain_length.min/max | exact |  | common/drc_rules.keywords:31; pcbnew/drc/drc_constraint_types.h:49 |
| 31 | `stub_length` | `NET_CHAIN_STUB_LENGTH_CONSTRAINT` | stub_length.max | exact |  | common/drc_rules.keywords:33; pcbnew/drc/drc_constraint_types.h:50 |
| 32 | `return_path` | `NET_CHAIN_RETURN_PATH_CONSTRAINT` | return_path.reference_layer, .reference_net | exact |  | common/drc_rules.keywords:32; common/drc_rules.keywords:41 |
| 33 | `skew` | `SKEW_CONSTRAINT` | skew.min/opt/max + skew.within_diff_pairs | exact |  | common/drc_rules.keywords:53; common/drc_rules.keywords:76 |
| 34 | `via_count` | `VIA_COUNT_CONSTRAINT` | via_count.min/max (count) | exact |  | common/drc_rules.keywords:72; pcbnew/drc/drc_constraint_types.h:56 |
| 35 | `diff_pair_gap` | `DIFF_PAIR_GAP_CONSTRAINT` | diff_pair_gap.min/opt/max | exact |  | common/drc_rules.keywords:14; pcbnew/drc/drc_constraint_types.h:53 |
| 36 | `diff_pair_uncoupled` | `MAX_UNCOUPLED_CONSTRAINT` | diff_pair_uncoupled.max | exact |  | common/drc_rules.keywords:15; pcbnew/drc/drc_constraint_types.h:54 |
| 37 | (no keyword) | `DIFF_PAIR_INTRA_SKEW_CONSTRAINT` | skew with within_diff_pairs=true | exact | API-only alias, no DRU keyword | pcbnew/drc/drc_constraint_types.h:55; pcbnew/drc/drc_rule.cpp:191-197 |
| 38 | `via_dangling` | `VIA_DANGLING_CONSTRAINT` | dangling_via (severity-only rule) | exact |  | common/drc_rules.keywords:73; pcbnew/drc/drc_constraint_types.h:62 |
| 39 | `microvia_stack_depth` | `MICROVIA_STACK_DEPTH_CONSTRAINT` | microvia_stack[depth].max_depth | exact |  | common/drc_rules.keywords:38; pcbnew/drc/drc_constraint_types.h:65 |
| 40 | `hole_aspect_ratio` | `HOLE_ASPECT_RATIO_CONSTRAINT` | hole_aspect_ratio.max (stored x 1000 in KiCad) | exact |  | common/drc_rules.keywords:23; pcbnew/drc/drc_constraint_types.h:66 |
| 41 | `assertion` | `ASSERTION_CONSTRAINT` | custom_assertion.expr; editor patterns recognised on import: permitted layers -> permitted_layers, allowed orientation -> component_orientation, enclosedByArea -> room_confinement, Height bounds -> component_height | exact | expression translated per section 2.4 | common/drc_rules.keywords:2; pcbnew/drc/drc_constraint_types.h:59 |
| 42 | `mechanical_clearance` | `PHYSICAL_CLEARANCE_CONSTRAINT` | physical_clearance.min | exact | deprecated spelling | common/drc_rules.keywords:35; pcbnew/drc/drc_rule_parser.cpp:520-524 |
| 43 | `mechanical_hole_clearance` | `PHYSICAL_HOLE_CLEARANCE_CONSTRAINT` | physical_hole_clearance.min | exact | deprecated spelling | common/drc_rules.keywords:36; pcbnew/drc/drc_rule_parser.cpp:525-529 |
| 44 | `hole` | `HOLE_SIZE_CONSTRAINT` | hole_size | exact | deprecated spelling | common/drc_rules.keywords:22; pcbnew/drc/drc_rule_parser.cpp:530-534 |
| 45 | (no keyword) | `NULL_CONSTRAINT` | none | none: internal 'no rule matched' result, not a rule |  | pcbnew/drc/drc_constraint_types.h:25; pcbnew/drc/drc_engine.cpp:2234-2238 |

#### 4.1.2 DRU clauses, evaluation semantics and implicit rules

Rows 19-40 are evaluation semantics; each points at the precedence rule of section 5 that reproduces it.

| # | Kind | Clause / topic | Canonical target | Quality | Note | Provenance |
|---|---|---|---|---|---|---|
| 1 | top-level clause | version | @ruleset.version | model | exporter writes the version KiCad expects | pcbnew/dialogs/panel_setup_rules_help_1clauses.md:3; pcbnew/drc/drc_rule_parser.cpp:187-216 |
| 2 | top-level clause | rule | @id + @native.name | model |  | pcbnew/dialogs/panel_setup_rules_help_1clauses.md:5; pcbnew/drc/drc_rule_parser.cpp:218-227 |
| 3 | top-level clause (separate grammar) | assign_component_class | @scope component_class() membership source | model | component-class assignment, not a DRC rule | common/drc_rules.keywords:3; pcbnew/drc/drc_rule_parser.cpp:238-300 |
| 4 | rule clause | constraint | @kind + @constraint | model |  | pcbnew/dialogs/panel_setup_rules_help_1clauses.md:11; pcbnew/drc/drc_rule_parser.cpp:330-332 |
| 5 | rule clause | condition | @scope.expr (+ expr_b) | model | section 2.4 | pcbnew/dialogs/panel_setup_rules_help_1clauses.md:13; pcbnew/drc/drc_rule_parser.cpp:334-360 |
| 6 | rule clause | layer | @layers | model | outer, inner and wildcards kept as layer tokens | pcbnew/dialogs/panel_setup_rules_help_1clauses.md:15; pcbnew/drc/drc_rule_parser.cpp:362-367 |
| 7 | rule clause | severity | @severity | model | section 2.7 | pcbnew/dialogs/panel_setup_rules_help_1clauses.md:17; pcbnew/drc/drc_rule_parser.cpp:369-371 |
| 8 | rule clause (error text) | disallow | none | none: parser error text only |  | pcbnew/drc/drc_rule_parser.cpp:377-378; pcbnew/drc/drc_rule_parser.cpp:328-379 |
| 9 | constraint argument | min / opt / max | @constraint.min/opt/max | model |  | pcbnew/drc/drc_rule_parser.cpp:809-870; pcbnew/drc/drc_rule_parser.cpp:483-513 |
| 10 | constraint argument | units | @constraint.unit | model | accepted units are mil, mm, in, deg, fs and ps: mil and in are converted to mm, deg kept, fs converted to ps; a value without a unit gets a 'Missing units' warning and is read as internal units (nm for lengths), except unitless constraints (counts, ratio) and `track_angle` (read as deg) | pcbnew/pcbexpr_evaluator.cpp:857-863; common/libeval_compiler/libeval_compiler.cpp:1093-1104; pcbnew/drc/drc_rule_parser.cpp:612-623 |
| 11 | constraint argument | within_diff_pairs | skew.within_diff_pairs | exact |  | common/drc_rules.keywords:76; pcbnew/drc/drc_rule_parser.cpp:798-807 |
| 12 | constraint argument | disallow kinds | disallow.kinds | exact | one row per token in the table below this one | pcbnew/drc/drc_rule_parser.cpp:630-671; pcbnew/drc/drc_constraint_types.h:70-98 |
| 13 | constraint argument | zone_connection values | zone_connection.value | exact |  | pcbnew/drc/drc_rule_parser.cpp:709-736 |
| 14 | constraint argument | return_path sub-clauses | return_path.reference_layer, .reference_net | exact |  | pcbnew/drc/drc_rule_parser.cpp:672-708 |
| 15 | constraint argument (help only) | time_domain | none | none: help-only token; time domain is chosen by ps or fs units |  | pcbnew/dialogs/panel_setup_rules_help_5examples.md:38; pcbnew/dialogs/panel_setup_rules_help_6notes.md:34 |
| 16 | lexical | comments | none | none: comment syntax |  | pcbnew/dialogs/panel_setup_rules_help_5examples.md:17; pcbnew/dialogs/panel_setup_rules_help_6notes.md:11 |
| 17 | lexical | text variables | @native.text (variables expanded on import) | model |  | pcbnew/drc/drc_engine.cpp:693-709; pcbnew/drc/drc_rule.cpp:32-35 |
| 18 | static check | B in unary constraints | arity check (expr_b forbidden for unary kinds) | model | section 2.2 | pcbnew/drc/drc_rule_parser.cpp:385-406; pcbnew/drc/drc_rule.cpp:88-118 |
| 19 | evaluation | load order | section 5, P2 | model |  | pcbnew/drc/drc_engine.cpp:906-927; pcbnew/drc/drc_engine.cpp:718-719 |
| 20 | evaluation | per-type lists | section 5, P2 | model |  | pcbnew/drc/drc_engine.cpp:744-767 |
| 21 | evaluation | matching | section 5, P3 | model |  | pcbnew/drc/drc_engine.cpp:1997-2003; pcbnew/drc/drc_engine.cpp:1077-1137 |
| 22 | evaluation | commutative conditions | section 5, P4 | model |  | pcbnew/drc/drc_rule_condition.cpp:105-117 |
| 23 | evaluation | layer test | section 5, P4 | model |  | pcbnew/drc/drc_engine.cpp:1819-1835 |
| 24 | evaluation | condition cache | none | none: evaluation cache, no semantic effect |  | pcbnew/drc/drc_engine.cpp:1413-1418; pcbnew/drc/drc_engine.cpp:1639-1665 |
| 25 | evaluation | net ties | section 5, P10 | model |  | pcbnew/drc/drc_engine.cpp:1142-1178 |
| 26 | evaluation | local overrides (clearance, hole_clearance) | section 5, P9 | model |  | pcbnew/drc/drc_engine.cpp:1180-1250 |
| 27 | evaluation | local overrides (thermal, mask, paste, zone connection) | section 5, P9 | model |  | pcbnew/drc/drc_engine.cpp:1251-1399 |
| 28 | evaluation | explicit beats implicit | section 5, P8 | model |  | pcbnew/drc/drc_engine.cpp:2025-2026 |
| 29 | evaluation | parent-footprint inheritance | section 5, P9 | model |  | pcbnew/drc/drc_engine.cpp:2028-2060 |
| 30 | evaluation | fallbacks after rules: board setup for mask and paste (:2062-2083, pads in footprints only); footprint local zone connection, then the zone's own connection, thermal gap and spoke width (:2164-2232) | section 5, P9 (origin item_default for the footprint and zone values) | model |  | pcbnew/drc/drc_engine.cpp:2062-2083; pcbnew/drc/drc_engine.cpp:2164-2232 |
| 31 | evaluation | clearance post-processing | section 5, P9 | model |  | pcbnew/drc/drc_engine.cpp:2086-2162 |
| 32 | evaluation | netclass clearance ordering | section 5, P6 | model |  | pcbnew/drc/drc_engine.cpp:456-468; pcbnew/drc/drc_engine.cpp:1670-1697 |
| 33 | evaluation | netclass fast path | section 5, P6 | model |  | pcbnew/drc/drc_engine.cpp:1929-1994 |
| 34 | evaluation | diff_pair_gap vs clearance | section 5, P7 | model |  | pcbnew/drc/drc_engine.cpp:2005-2022 |
| 35 | evaluation | clearance cap | section 5, P11 | model |  | pcbnew/drc/drc_engine.cpp:1110-1120 |
| 36 | evaluation | disallow matching | section 5, P12 | model |  | pcbnew/drc/drc_engine.cpp:1716-1817 |
| 37 | evaluation | assertions | section 5, P12 | model |  | pcbnew/drc/drc_engine.cpp:2278-2347; pcbnew/drc/drc_engine.cpp:2309-2313 |
| 38 | evaluation | length/skew grouping | section 5, P12 | model |  | pcbnew/drc/drc_test_provider_matched_length.cpp:393-405 |
| 39 | evaluation | teardrops | section 2.3 (teardrops count as tracks) | model |  | pcbnew/drc/drc_rule_condition.cpp:98-103; pcbnew/drc/drc_engine.cpp:1762-1767 |
| 40 | evaluation | no match | section 5, P3 | model |  | pcbnew/drc/drc_engine.cpp:2234-2238 |
| 41 | implicit rule | board setup constraints | track_width, connection_width, hole_size, annular_width, via_diameter, hole_to_hole [origin board_floor] | exact |  | pcbnew/drc/drc_engine.cpp:170-195; pcbnew/board_design_settings.cpp:296-326 |
| 42 | implicit rule | board setup constraints zone fill strategy | starved_thermal.min_resolved [board_floor] | exact |  | pcbnew/drc/drc_engine.cpp:197-201; pcbnew/board_design_settings.cpp:344 |
| 43 | implicit rule | board setup constraints silk / silk text height / silk text thickness | silk_clearance, text_height, text_thickness [board_floor, silk layers] | exact |  | pcbnew/drc/drc_engine.cpp:203-221; pcbnew/board_design_settings.cpp:332-353 |
| 44 | implicit rule | board setup constraints hole / edge / courtyard | hole_clearance, edge_clearance, courtyard_clearance [board_floor] | lossy: courtyard value unset in source (Q-K5) |  | pcbnew/drc/drc_engine.cpp:223-236; pcbnew/board_design_settings.cpp:328-359 |
| 45 | implicit rule | board setup constraints micro-via | hole_size, via_diameter [board_floor, scope type(via.micro)] | exact |  | pcbnew/drc/drc_engine.cpp:240-251; pcbnew/board_design_settings.cpp:316-322 |
| 46 | implicit rule | barcode visual separation default | physical_clearance.min = 1.0 mm [tool_default, scope type(barcode)] | exact |  | pcbnew/drc/drc_engine.cpp:253-260 |
| 47 | implicit rule | netclass '&lt;name>' | clearance, track_width, diff_pair_gap, via_diameter, hole_size [origin netclass] | exact |  | pcbnew/drc/drc_engine.cpp:262-454 |
| 48 | implicit rule | netclass '&lt;name>' diff pair (gap < clearance) | clearance [netclass, pair coupled_diffpair] | exact |  | pcbnew/drc/drc_engine.cpp:353-371 |
| 49 | implicit rule | tuning profile '&lt;name>' | track_width, diff_pair_gap, clearance per layer [netclass, from tuning profile] + impedance_target | lossy: ohm target not kept |  | pcbnew/drc/drc_engine.cpp:473-598 |
| 50 | implicit rule | keepout area ['&lt;name>'] [of &lt;ref>] | disallow [explicit, scope in_area(uuid, intersects_fill)] | exact |  | pcbnew/drc/drc_engine.cpp:600-678 |
| 51 | constraint options | `DRC_CONSTRAINT::OPTIONS`: `SKEW_WITHIN_DIFF_PAIRS`, `SPACE_DOMAIN`, `TIME_DOMAIN` | SKEW_WITHIN_DIFF_PAIRS -> skew.within_diff_pairs (@constraint.options); SPACE_DOMAIN / TIME_DOMAIN -> @constraint.unit mm or ps | model | the parser sets the domain from the units of the values (time units select TIME_DOMAIN) and sets SKEW_WITHIN_DIFF_PAIRS from `(within_diff_pairs)`; the enum also has NUM_OPTIONS, a size sentinel | pcbnew/drc/drc_rule.h:112-120; pcbnew/drc/drc_rule_parser.cpp:500-512, :798-800 |
| 52 | implicit-rule source | `DRC_IMPLICIT_SOURCE`: `NONE`, `BOARD_SETUP_CONSTRAINT`, `BARCODE_DEFAULTS`, `KEEPOUT`, `NET_CLASS`, `TUNING_PROFILE` | @origin: NONE (user DRU rule) -> explicit; BOARD_SETUP_CONSTRAINT -> board_floor; BARCODE_DEFAULTS -> tool_default; KEEPOUT -> explicit; NET_CLASS -> netclass; TUNING_PROFILE -> netclass (from tuning profile) | model | the engine uses NET_CLASS to order netclass clearance rules (:1672, :1895) | pcbnew/drc/drc_constraint_types.h:86-94; pcbnew/drc/drc_rule.cpp:45; pcbnew/drc/drc_engine.cpp:171-241, :256, :283-433, :484-549, :612-631 |

`disallow` argument tokens (row 12). The parser ORs one `DRC_DISALLOW_*` flag per token into the constraint
(pcbnew/drc/drc_rule_parser.cpp:639-653; flags pcbnew/drc/drc_constraint_types.h:72-82):

| # | DRU token | `DRC_DISALLOW_*` flag(s) | Canonical `disallow.kinds` member(s) | Exporter writes | Quality | Provenance |
|---|---|---|---|---|---|---|
| D1 | `track` | `DRC_DISALLOW_TRACKS` | track | `track` | exact | pcbnew/drc/drc_rule_parser.cpp:639; pcbnew/drc/drc_constraint_types.h:76; common/drc_rules.keywords:66 |
| D2 | `via` | `DRC_DISALLOW_THROUGH_VIAS`, `DRC_DISALLOW_BLIND_VIAS`, `DRC_DISALLOW_BURIED_VIAS`, `DRC_DISALLOW_MICRO_VIAS` | through_via, blind_via, buried_via, micro_via (the canonical `via` member is shorthand for these four) | `via` when all four via members are set; otherwise one token per member (rows D3-D6) | exact | pcbnew/drc/drc_rule_parser.cpp:640-643; pcbnew/drc/drc_constraint_types.h:72-75; common/drc_rules.keywords:71 |
| D3 | `through_via` | `DRC_DISALLOW_THROUGH_VIAS` | through_via | `through_via` (KiCad 10.0 and later; for 9.x see section 4.1.10 row 16) | exact | pcbnew/drc/drc_rule_parser.cpp:644; pcbnew/drc/drc_constraint_types.h:72; common/drc_rules.keywords:65 |
| D4 | `blind_via` | `DRC_DISALLOW_BLIND_VIAS` | blind_via | `blind_via` (KiCad 10.0 and later; section 4.1.10 row 10) | exact | pcbnew/drc/drc_rule_parser.cpp:645; pcbnew/drc/drc_constraint_types.h:74; common/drc_rules.keywords:6 |
| D5 | `buried_via` | `DRC_DISALLOW_BURIED_VIAS` | buried_via | `buried_via` | exact | pcbnew/drc/drc_rule_parser.cpp:646; pcbnew/drc/drc_constraint_types.h:75; common/drc_rules.keywords:7 |
| D6 | `micro_via` | `DRC_DISALLOW_MICRO_VIAS` | micro_via | `micro_via` | exact | pcbnew/drc/drc_rule_parser.cpp:647; pcbnew/drc/drc_constraint_types.h:73; common/drc_rules.keywords:37 |
| D7 | `pad` | `DRC_DISALLOW_PADS` | pad | `pad` | exact | pcbnew/drc/drc_rule_parser.cpp:648; pcbnew/drc/drc_constraint_types.h:77; common/drc_rules.keywords:46 |
| D8 | `zone` | `DRC_DISALLOW_ZONES` | zone_fill (renamed: the canonical member names what is forbidden, copper fill, as KiCad's `Keep_Out_Zone_Fills` property does) | `zone` | exact | pcbnew/drc/drc_rule_parser.cpp:649; pcbnew/drc/drc_constraint_types.h:78; common/drc_rules.keywords:77 |
| D9 | `text` | `DRC_DISALLOW_TEXTS` | text | `text` | exact | pcbnew/drc/drc_rule_parser.cpp:650; pcbnew/drc/drc_constraint_types.h:79; common/drc_rules.keywords:59 |
| D10 | `graphic` | `DRC_DISALLOW_GRAPHICS` | graphic | `graphic` | exact | pcbnew/drc/drc_rule_parser.cpp:651; pcbnew/drc/drc_constraint_types.h:80; common/drc_rules.keywords:21 |
| D11 | `hole` | `DRC_DISALLOW_HOLES` | hole | `hole` | exact | pcbnew/drc/drc_rule_parser.cpp:652; pcbnew/drc/drc_constraint_types.h:81; common/drc_rules.keywords:22 |
| D12 | `footprint` | `DRC_DISALLOW_FOOTPRINTS` | footprint | `footprint` | exact | pcbnew/drc/drc_rule_parser.cpp:653; pcbnew/drc/drc_constraint_types.h:82; common/drc_rules.keywords:20 |

The header also defines two composite masks, `DRC_DISALLOW_VIAS` (= through vias only) and `DRC_DISALLOW_BB_VIAS`
(blind or buried), which the parser does not produce (pcbnew/drc/drc_constraint_types.h:97-98). A string argument
(`"via"`) is accepted as the same token (pcbnew/drc/drc_rule_parser.cpp:634-635).

#### 4.1.3 Board-setup keys in `.kicad_pro` (`board.design_settings`)

| # | JSON key | Unit | Default | Canonical target | Quality | Note | Provenance |
|---|---|---|---|---|---|---|---|
| 1 | `rules.min_clearance` | mm | 0.0 | clearance.min [board_floor] | exact | KiCad floors implicit and local results only (section 5, P8) | pcbnew/board_design_settings.cpp:152, :292-294; include/board_design_settings_defaults.h:69 |
| 2 | `rules.min_connection` | mm | 0.0 | connection_width.min [board_floor] | exact |  | pcbnew/board_design_settings.cpp:153, :296-298; include/board_design_settings_defaults.h:70 |
| 3 | `rules.min_track_width` | mm | 0.2 | track_width.min [board_floor] | exact |  | pcbnew/board_design_settings.cpp:154, :300-302; include/board_design_settings_defaults.h:71 |
| 4 | `rules.min_via_annular_width` | mm | 0.1 from the constructor, (0.5 − 0.3) / 2. The `PARAM` default argument is 0.5 (Q-K4) | annular_width.min [board_floor] | exact | two defaults in source (Q-K4) | pcbnew/board_design_settings.cpp:155, :304-306; pcbnew/drc/drc_engine.cpp:185-187 |
| 5 | `rules.min_via_diameter` | mm | 0.5 | via_diameter.min [board_floor] | exact | advisory: ineffective below drill + 2 x annulus | pcbnew/board_design_settings.cpp:156, :308-310, :1761-1773; include/board_design_settings_defaults.h:72 |
| 6 | `rules.min_through_hole_diameter` | mm | 0.3 | hole_size.min [board_floor] | exact |  | pcbnew/board_design_settings.cpp:157, :312-314; include/board_design_settings_defaults.h:73 |
| 7 | `rules.min_microvia_diameter` | mm | 0.2 | via_diameter.min [board_floor, type(via.micro)] | exact |  | pcbnew/board_design_settings.cpp:158, :316-318; include/board_design_settings_defaults.h:74 |
| 8 | `rules.min_microvia_drill` | mm | 0.1 | hole_size.min [board_floor, type(via.micro)] | exact |  | pcbnew/board_design_settings.cpp:159, :320-322; include/board_design_settings_defaults.h:75 |
| 9 | `rules.min_hole_to_hole` | mm | 0.25 | hole_to_hole.min [board_floor] | exact |  | pcbnew/board_design_settings.cpp:162, :324-326; include/board_design_settings_defaults.h:76 |
| 10 | `rules.min_hole_clearance` | mm | 0.25 | hole_clearance.min [board_floor] | exact |  | pcbnew/board_design_settings.cpp:161, :328-330; include/board_design_settings_defaults.h:77 |
| 11 | `rules.min_silk_clearance` | mm | 0.0 | silk_clearance.min + silk_to_mask_clearance.min + silk_to_edge_clearance.min [board_floor] | exact | one KiCad value feeds three canonical kinds; export takes the smallest of the three | pcbnew/board_design_settings.cpp:163, :332-334; include/board_design_settings_defaults.h:82 |
| 12 | `rules.min_groove_width` | mm | 0.0 | creepage.min_groove_width | exact |  | pcbnew/board_design_settings.cpp:167, :336-338; include/board_design_settings.h:661 |
| 13 | `rules.min_resolved_spokes` | count | 2 | starved_thermal.min_resolved [board_floor] | exact |  | pcbnew/board_design_settings.cpp:164, :340-345; include/board_design_settings_defaults.h:85 |
| 14 | `rules.min_text_height` | mm | 0.8 (1.0 × 0.8) | text_height.min [board_floor, silk layers] | exact |  | pcbnew/board_design_settings.cpp:165, :347-349; include/board_design_settings_defaults.h:31 |
| 15 | `rules.min_text_thickness` | mm | 0.08 (0.1 × 0.8) | text_thickness.min [board_floor, silk layers] | exact |  | pcbnew/board_design_settings.cpp:166, :351-353; include/board_design_settings_defaults.h:35 |
| 16 | `rules.min_copper_edge_clearance` | mm | 0.5 | edge_clearance.min [board_floor] | lossy: the -0.01 mm legacy flag has no canonical value |  | pcbnew/board_design_settings.cpp:160, :355-359; include/board_design_settings_defaults.h:79-81 |
| 17 | `rules.max_error` | mm | `ARC_HIGH_DEF` (constant not in checkout; every golden project stores 0.005) | drc_engine_settings.arc_tolerance | exact |  | pcbnew/board_design_settings.cpp:221, :1136-1138, :1682-1684; include/board_design_settings_defaults.h:87-88 |
| 18 | `rules.solder_mask_to_copper_clearance` | mm | 0.0 | solder_mask_bridge.min_copper_to_aperture | exact |  | pcbnew/board_design_settings.cpp:232, :1140-1142; include/board_design_settings.h:693-694 |
| 19 | `rules.use_height_for_length_calcs` | bool | true | length.include_via_height | lossy: consumer not in checkout (Q-K8) |  | pcbnew/board_design_settings.cpp:223, :289-290 |
| 20 | `rule_severities` | map: settings key → `error` / `warning` / `ignore` | see section 7 | @severity per kind (rule-set default) via section 4.1.12 | model |  | pcbnew/board_design_settings.cpp:169-219, :361-400 |
| 21 | `drc_exclusions` | array of exclusion objects (legacy strings migrated in schema 2→3) | `[]` | @waiver (finding status dismissed) | model |  | pcbnew/board_design_settings.cpp:402-430, :1483-1521 |
| 22 | `track_widths` | mm list | `[]`; golden projects store 0.0 first, the netclass placeholder | router_presets.widths | exact |  | pcbnew/board_design_settings.cpp:432-457; include/board_design_settings.h:637 |
| 23 | `via_dimensions` | list of {`diameter`, `drill`} in mm | `[]` | router_presets.vias | exact |  | pcbnew/board_design_settings.cpp:459-497, :1723-1742 |
| 24 | `diff_pair_dimensions` | list of {`width`, `gap`, `via_gap`} in mm | `[]` | router_presets.pairs | exact |  | pcbnew/board_design_settings.cpp:609-653, :1712-1721 |
| 25 | `via_stack_presets` | {`name`, `start_layer`, `end_layer`, `staggered`, `via_size` mm, `via_drill` mm, `use_netclass`, `filled`, `capped`, `pitch` mm} | `[]` | router_presets.via_stacks | exact |  | pcbnew/board_design_settings.cpp:499-607; include/board_design_settings_defaults.h:42-43 |
| 26 | `teardrop_options, teardrop_parameters` | target flags (bool); `td_maxlen`, `td_maxheight` mm; ratios | none set in BDS | teardrop_style | exact |  | pcbnew/board_design_settings.cpp:656-798, :1697-1710 |
| 27 | `tuning_pattern_settings` | `min_amplitude`, `max_amplitude`, `spacing` mm; `corner_radius_percentage` %; `single_sided` bool | spacing 0.6 mm single and skew, 1.0 mm diff pair | meander_style | exact |  | pcbnew/board_design_settings.cpp:259-261, :800-866; include/board_design_settings_defaults.h:66-67 |
| 28 | `defaults.zones.*` | `min_clearance`, `min_thickness`, `thermal_relief_gap`, `thermal_relief_spoke_width`, `hatch_*`, `border_hatch_pitch` mm; `min_island_area` mm²; `pad_connection`, `remove_islands`, `fill_mode`, `corner_smoothing` enums | `ZONE_*_MM` constants (not in checkout); `pad_connection` THERMAL; `remove_islands` ALWAYS; `min_island_area` 10.0 mm² | zone defaults [tool_default]: clearance, zone_min_width, thermal_relief_gap, thermal_spoke_width, zone_connection, zone_fill_style | exact |  | pcbnew/board_design_settings.cpp:1010-1134 |
| 29 | `zones_allow_external_fillets` | bool | false (forced true for projects with a `legacy` section) | zone_fill_style.external_fillets | exact |  | pcbnew/board_design_settings.cpp:222, :1144-1145, :1580-1584 |
| 30 | `defaults.*_line_width, defaults.*_text_size_v/_h, defaults.*_text_thickness` | mm | silk line 0.1; copper line 0.2; edge 0.05; courtyard 0.05; fab and other 0.1. Text size / stroke: silk 1.0 / 0.1; copper 1.5 / 0.3; others 1.0 / 0.15 | none | none: creation defaults for new objects, not constraints |  | include/board_design_settings_defaults.h:25-37, :93-94; pcbnew/board_design_settings.cpp:83-125, :873-968 |

#### 4.1.4 Board-file `(setup)` values

| # | Member | Unit | Default | Canonical target | Quality | Note | Provenance |
|---|---|---|---|---|---|---|---|
| 1 | `m_SolderMaskExpansion` | mm | 0.0 | solder_mask_expansion.opt [board default] | exact | token pad_to_mask_clearance | include/board_design_settings_defaults.h:50-52; pcbnew/board_design_settings.cpp:230, :1156-1157, :1686-1687; token: pcbnew/pcb_io/kicad_sexpr/pcb_io_kicad_sexpr_parser.cpp:3297-3300† |
| 2 | `m_SolderMaskMinWidth` | mm | 0.0 | solder_mask_bridge.min_web [board default] + solder_mask_sliver.min [board_floor] | lossy: KiCad uses the value only as the bridge criterion between different-net apertures (apertures inflated by half the web; aperture-to-aperture collision distance; same-net pairs filtered out first) and never reports a sliver; the canonical import keeps it as the web of solder_mask_bridge, which is what KiCad enforces, and as the sliver floor, which the kernel checks more strictly than KiCad; export writes solder_mask_sliver.min | token solder_mask_min_width (a missing token means 0, sexpr parser :3042-3044†) | include/board_design_settings.h:691-692; include/board_design_settings_defaults.h:56; pcbnew/drc/drc_test_provider_solder_mask.cpp:77, :205, :506, :584-587, :840, :1003; token: pcbnew/pcb_io/kicad_sexpr/pcb_io_kicad_sexpr_parser.cpp:3302-3305† |
| 3 | `m_SolderPasteMargin` | mm | 0.0 | solder_paste_margin.abs [board default] | exact |  | include/board_design_settings.h:696-698; include/board_design_settings_defaults.h:58 |
| 4 | `m_SolderPasteMarginRatio` | fraction (DRC holds it × 1000) | 0.0 | solder_paste_margin.ratio [board default] | exact |  | include/board_design_settings_defaults.h:59; pcbnew/board_design_settings.cpp:239, :1165-1166, :1694-1695 |
| 5 | `m_AllowSoldermaskBridgesInFPs` | bool | false | solder_mask_bridge.allow_in_footprint | exact |  | include/board_design_settings.h:699; pcbnew/board_design_settings.cpp:241 |
| 6 | `m_TentViasFront, m_TentViasBack` | bool | true, true | via_protection.tent_top/tent_bottom | exact |  | include/board_design_settings.h:701-702; pcbnew/board_design_settings.cpp:242-243 |
| 7 | `m_CoverViasFront/Back, m_PlugViasFront/Back, m_CapVias` | bool | false | via_protection.cover_*/plug_*/cap | exact |  | include/board_design_settings.h:704-711; pcbnew/board_design_settings.cpp:245-251 |
| 8 | `m_FillVias` | bool | false | via_protection.fill (+ microvia_stack[not_filled] fallback) | exact |  | include/board_design_settings.h:713-714; pcbnew/board_design_settings.cpp:253 |
| 9 | `m_boardThickness` | mm | 1.6 | stackup.board_thickness | exact |  | include/board_design_settings_defaults.h:45-46; include/board_design_settings.h:782 |

#### 4.1.5 Hard-coded implicit rules and item-level overrides

| # | Behaviour | Canonical target | Quality | Note | Provenance |
|---|---|---|---|---|---|
| 1 | Implicit rules built from board setup, in this order: "board setup constraints", "… zone fill strategy", "… silk", "… silk text height", "… silk text thickness", "… hole", "… edge", "… courtyard", "… micro-via" | track_width, connection_width, hole_size, annular_width, via_diameter, hole_to_hole, starved_thermal, silk_clearance, text_height, text_thickness, hole_clearance, edge_clearance, courtyard_clearance [board_floor] | exact |  | pcbnew/drc/drc_engine.cpp:163-251 |
| 2 | Barcode visual separation, condition `A.Type == 'Barcode'` | physical_clearance.min = 1.0 mm [tool_default, type(barcode)] | exact |  | pcbnew/drc/drc_engine.cpp:253-260 |
| 3 | Courtyard implicit rule. Line 235 sets the min on `holeToHoleConstraint`, not on the courtyard constraint, so the courtyard value is left unset | courtyard_clearance [tool_default] | lossy: value unset in source (Q-K5) |  | pcbnew/drc/drc_engine.cpp:233-236 |
| 4 | Per-netclass rules, condition `A.hasExactNetclass('<name>')`, plus `A.inDiffPair('*')`, `A.Via_Type == / != 'Micro'`, `AB.isCoupledDiffPair()` variants | clearance, track_width, diff_pair_gap, via_diameter, hole_size [origin netclass] | exact |  | pcbnew/drc/drc_engine.cpp:262-471 |
| 5 | Netclass clearance rules are sorted by min ascending. Selection is "winner takes all", so local clearances are max-ed in code, not by rules | section 5, P6 (larger netclass clearance wins) | model |  | pcbnew/drc/drc_engine.cpp:456-468, :2086-2088 |
| 6 | A netclass whose name cannot be quoted as a string literal (a name ending in a backslash) gets no implicit rules; the misc provider reports it | none | none: KiCad quoting limitation; canonical strings escape `'` as `\'` and can hold such a name, but it cannot be exported to KiCad (section 2.4) |  | pcbnew/drc/drc_engine.cpp:270-274; common/libeval_compiler/libeval_compiler.cpp:1472-1480; pcbnew/drc/drc_test_provider_misc.cpp:621-645 |
| 7 | Tuning-profile rules: per layer, TRACK_WIDTH with min = opt = max = profile width (min floored at `rules.min_track_width`); for differential profiles also DIFF_PAIR_GAP (floored at `rules.min_clearance`) and a CLEARANCE override when the gap is below the netclass clearance. Severity comes from key `tuning_profile_track_geometries` (default ignore) | track_width, diff_pair_gap, clearance [netclass, from tuning profile]; severity from tuning_profile_consistency[implicit_geometry] | exact |  | pcbnew/drc/drc_engine.cpp:473-598; pcbnew/board_design_settings.cpp:217 |
| 8 | Keepout / rule areas become implicit rules, condition `A.intersectsKeepout('<uuid>')`, layer set = zone layers, flags tracks / vias / pads / zones / footprints | disallow [explicit, scope in_area(uuid, intersects_fill), layers = zone layers] | exact |  | pcbnew/drc/drc_engine.cpp:600-660 |
| 9 | Footprint net ties: for two copper items where one is connected and one is not, or both are in the same footprint, CLEARANCE becomes 0 mm when the connected item's net is in the footprint's net-tie set for the other item | net_tie | exact |  | pcbnew/drc/drc_engine.cpp:1142-1176 |
| 10 | Item clearance overrides (pad, or the pad's parent footprint) win over rules but are floored at `rules.min_clearance` (copper) or `rules.min_hole_clearance` (holes). Only PAD implements the override; other connected items, zones included, return none, so a zone's clearance takes the row 11 path | clearance, hole_clearance [item_override], floored at board_floor | exact | section 5, P9 | pcbnew/drc/drc_engine.cpp:1180-1245; pcbnew/pad.cpp:1913-1921†; pcbnew/board_connected_item.h:169-172†; pcbnew/zone.h:182, :191† |
| 11 | Item local clearance is max-ed with the rule result, then floored at `rules.min_clearance` for copper pairs | clearance [item_override on zones: max with rule] | exact | section 5, P9 | pcbnew/drc/drc_engine.cpp:2089-2160 |
| 12 | Pad zone-connection override (unless INHERITED), checked before rules; the footprint's own zone connection is checked only after rules on the pad and on the footprint, then the zone's pad-connection setting | zone_connection: pad value [item_override]; footprint value and zone value [item_default, in that order] | exact | section 5, P9 | pcbnew/drc/drc_engine.cpp:1251-1266, :2028-2060, :2164-2197 |
| 13 | Pad thermal-gap override (when > 0), checked before rules; after rules on the pad and on the footprint, the zone's own thermal gap | thermal_relief_gap: pad value [item_override]; zone value [item_default] | exact | section 5, P9 | pcbnew/drc/drc_engine.cpp:1268-1283, :2199-2213 |
| 14 | Pad spoke-width override (when > 0), raised to the zone min thickness, checked before rules; after rules on the pad and on the footprint, the zone's own spoke width | thermal_spoke_width: pad value [item_override]; zone value [item_default] | exact | raised to the zone minimum width; section 5, P9 | pcbnew/drc/drc_engine.cpp:1285-1300, :2216-2230 |
| 15 | Local mask margin: pad, shape or track first, then the parent footprint's margin (both before rules), then rules on the item, rules on the parent footprint, then board setup (pads in footprints only) | solder_mask_expansion [item_override]; board value [board default] | exact | section 5, P9 | pcbnew/drc/drc_engine.cpp:1312-1340, :2028-2068 |
| 16 | Local paste margin and ratio: pad first, then footprint (both before rules), then rules on the item, rules on the parent footprint, then board setup | solder_paste_margin [item_override]; board value [board default] | exact | section 5, P9 | pcbnew/drc/drc_engine.cpp:1344-1399, :2028-2083 |
| 17 | DRC epsilon (tolerance) from advanced config `m_DRCEpsilon` (mm) | drc_engine_settings.epsilon | exact |  | pcbnew/board_design_settings.cpp:2054-2057; pcbnew/drc/drc_test_provider_copper_clearance.cpp:128 |
| 18 | Hole plating thickness from advanced config `m_HoleWallThickness` (mm) | drc_engine_settings.hole_wall_thickness | exact |  | pcbnew/board_design_settings.cpp:2060-2063 |

#### 4.1.6 Netclass fields (`net_settings.classes[]`)

| # | Netclass key | Unit | Default | Canonical target | Quality | Note | Provenance |
|---|---|---|---|---|---|---|---|
| 1 | `name` | string | `Default` | @scope netclass('&lt;name>') | model |  | common/netclass.cpp:31; common/project/net_settings.cpp:76, :131-133 |
| 2 | `priority` | integer; lower number wins | Default: 2147483647 (INT_MAX). Other classes: −1 when created. Label-only (undefined) classes: INT_MAX − 1. Old projects are migrated to 0, 1, 2… in list order | @netclass resolution (composite classes, per field) | model | section 5, P5 | common/netclass.cpp:57; common/project/net_settings.cpp:69, :77, :135-136, :607-616, :1025, :1195-1213 |
| 3 | `clearance` | mm | 0.2 | clearance.min [netclass] | exact |  | common/netclass.cpp:36, :66; common/project/net_settings.cpp:98-99, :141 |
| 4 | `track_width` | mm | 0.2 | track_width.opt [netclass] | exact | min stays the board floor | common/netclass.cpp:41, :69; common/project/net_settings.cpp:101-102, :144 |
| 5 | `via_diameter` | mm | 0.6 | via_diameter.opt [netclass, type(via) and not micro] | exact |  | common/netclass.cpp:37, :70; common/project/net_settings.cpp:104-105, :147 |
| 6 | `via_drill` | mm | 0.3 | hole_size.opt [netclass, type(via) and not micro] | exact |  | common/netclass.cpp:38, :67; common/project/net_settings.cpp:107-108, :150 |
| 7 | `microvia_diameter` | mm | 0.3 | via_diameter.opt [netclass, type(via.micro)] | exact |  | common/netclass.cpp:39, :71; common/project/net_settings.cpp:110-111, :153 |
| 8 | `microvia_drill` | mm | 0.1 | hole_size.opt [netclass, type(via.micro)] | exact |  | common/netclass.cpp:40, :68; common/project/net_settings.cpp:113-114, :156 |
| 9 | `diff_pair_width` | mm | 0.2 | track_width.opt [netclass, diffpair('*')] | exact |  | common/netclass.cpp:42, :72; common/project/net_settings.cpp:116-117, :159 |
| 10 | `diff_pair_gap` | mm | 0.25 | diff_pair_gap.opt [netclass] (+ clearance for coupled pairs when gap < clearance) | exact |  | common/netclass.cpp:43, :73; common/project/net_settings.cpp:119-120, :162 |
| 11 | `diff_pair_via_gap` | mm | 0.25 | diff_pair_gap.via_gap [netclass] | exact | router value; no DRC rule | common/netclass.cpp:44, :74; common/project/net_settings.cpp:122-123, :165 |
| 12 | `tuning_profile` | string | empty | impedance_target.profile + tuning_profile_consistency | lossy: ohm values live in the project tuning profiles |  | common/netclass.cpp:58; common/project/net_settings.cpp:80, :138-139, :646-651 |
| 13 | `wire_width` | mil | 6 | none | none: schematic display only |  | common/netclass.cpp:46, :76; common/project/net_settings.cpp:88-90, :168 |
| 14 | `bus_width` | mil | 12 | none | none: schematic display only |  | common/netclass.cpp:47, :77; common/project/net_settings.cpp:92-93, :171 |
| 15 | `line_style` | enum int (0 = solid) | 0 | none | none: schematic display only |  | common/netclass.cpp:49, :78; common/project/net_settings.cpp:95-96, :174-175 |
| 16 | `pcb_color` | RGBA string | unspecified | none | none: display only |  | common/netclass.cpp:61; common/project/net_settings.cpp:79, :177-178 |
| 17 | `schematic_color` | RGBA string | unspecified | none | none: display only |  | common/netclass.cpp:62; common/project/net_settings.cpp:78, :180-183 |
| 18 | `(description)` | string | "This is the default net class." | none | none: not saved |  | common/netclass.cpp:134-136; common/project/net_settings.cpp:68 |

#### 4.1.7 Netclass assignment and resolution

| # | Item | Canonical | Quality | Note | Provenance |
|---|---|---|---|---|---|
| 1 | `netclass_patterns` | @netclass membership (patterns) | model | matcher syntax unverified (Q-K7) | common/project/net_settings.cpp:364-408, :762-774, :1068-1087; golden-rules/bms-c1.txt:14 |
| 2 | `netclass_assignments` | @netclass membership (explicit assignments) | model |  | common/project/net_settings.cpp:323-359, :1048-1066 |
| 3 | `net_chain_classes`, `net_chain_netclasses` | net_chain_class() membership | model |  | common/project/net_settings.cpp:255-318, :1089-1090 |
| 4 | `net_colors` | none | none: display only |  | common/project/net_settings.cpp:224-250 |
| 5 | No match | @netclass default = 'Default' | model |  | common/project/net_settings.cpp:1035-1037, :1092-1139 |
| 6 | Unknown class name | @netclass implicit empty class | model |  | common/project/net_settings.cpp:1014-1033 |
| 7 | Several matches | @netclass composite resolution | model | section 5, P5 | common/project/net_settings.cpp:1141-1165, :1192-1240 |
| 8 | Composite classes in DRC | netclass-origin rules for composite classes | model |  | pcbnew/drc/drc_engine.cpp:453-454 |
| 9 | KiCad 5 board `(net_class <name> "<description>" ...)` | @netclass definitions (KiCad 5 keys) | model | legacy parser not in checkout (GK-K) | golden-rules/eez-psu-h24005.txt:37-38; golden-rules/esp32-evb.txt:102 |

#### 4.1.8 Where KiCad keeps rule data

| # | Store | Canonical | Quality | Note | Provenance |
|---|---|---|---|---|---|
| 1 | Project file `.kicad_pro` | @ruleset (board_floor, severities, presets) | model |  | pcbnew/board_design_settings.cpp:44-48; golden-rules/bms-c1.txt:3 (extracted from this path by scripts/golden_rules.py:23 |
| 2 | Project file `.kicad_pro` | @netclass definitions and membership | model |  | common/project/net_settings.cpp:41, :64-65, :189, :224, :255, :289, :323, :364; golden-rules/bms-c1.txt:11 (extracted by scripts/golden_rules.py:30) |
| 3 | Project file `.kicad_pro` | ERC rule set (severities, waivers, pin matrix) | model |  | eeschema/erc/erc_settings.cpp:32, :89-90, :137, :173, :203; golden-rules/bms-c1.txt:15-16 (extracted by scripts/golden_rules.py:34) |
| 4 | Board file `.kicad_pcb` | @ruleset (board defaults for mask, paste, via protection, stackup) | model |  | pcbnew/board_design_settings.cpp:1153-1157 (comment: mask and paste margins moved back to the board); golden-rules/bms-c1.txt:18 |
| 5 | Board file (KiCad 5 and earlier) | @netclass definitions (KiCad 5) | model |  | pcbnew/board_design_settings.cpp:54-58; golden-rules/eez-psu-h24005.txt:37 |

#### 4.1.9 DRU severity tokens

| # | DRU severity | Canonical | Provenance |
|---|---|---|---|
| 1 | `warning` | @severity=warning | pcbnew/dialogs/panel_setup_rules_help_4severity_names.md:3 |
| 2 | `error` | @severity=error | pcbnew/dialogs/panel_setup_rules_help_4severity_names.md:4 |
| 3 | `exclusion` | @tags+=waive (findings emitted as dismissed) | pcbnew/dialogs/panel_setup_rules_help_4severity_names.md:5 |
| 4 | `ignore` | @severity=ignore | pcbnew/dialogs/panel_setup_rules_help_4severity_names.md:6 |
| 5 | (no severity clause) | @severity=inherit | pcbnew/drc/drc_rule.cpp:44 |
| 6 | (resolution order) | section 2.6 (severity resolution) | pcbnew/pcb_marker.cpp:413-430† |
| 7 | (implicit tuning-profile rules) | tuning_profile_consistency[implicit_geometry] severity | pcbnew/drc/drc_engine.cpp:482 |

#### 4.1.10 Version gates (exporter targets)

The canonical importer accepts the union; the exporter must check the target KiCad version.

| # | Topic | Detail (abridged) | Effect on the mapping | Provenance |
|---|---|---|---|---|
| 1 | source pin | All line numbers refer to this commit. Files marked with a dagger were fetched from GitLab at the same commit because the sparse checkout lacks them. | recorded for the exporter (no mapping change) | tools/generate_drc_rules_vocabulary.py:41-44 |
| 2 | vocabulary method | The committed header lists, per release series, the tokens, functions ... | recorded for the exporter (no mapping change) | include/drc_rules_release_vocabulary.h:30 |
| 3 | token added after 10.0 | `hole_aspect_ratio` | export to KiCad 10.0 or older: not expressible (lossy); import: no effect | common/drc_rules.keywords:23 |
| 4 | token added after 10.0 | `microvia_stack_depth` | export to KiCad 10.0 or older: not expressible (lossy); import: no effect | common/drc_rules.keywords:38 |
| 5 | token added after 10.0 | `net` | export to KiCad 10.0 or older: not expressible (lossy); import: no effect | common/drc_rules.keywords:41 |
| 6 | token added after 10.0 | `net_chain_length` | export to KiCad 10.0 or older: not expressible (lossy); import: no effect | common/drc_rules.keywords:31 |
| 7 | token added after 10.0 | `return_path` | export to KiCad 10.0 or older: not expressible (lossy); import: no effect | common/drc_rules.keywords:32 |
| 8 | token added after 10.0 | `stub_length` | export to KiCad 10.0 or older: not expressible (lossy); import: no effect | common/drc_rules.keywords:33 |
| 9 | token added in 10.0 | `assign_component_class` | export to KiCad 9.x: not expressible (lossy) | include/drc_rules_release_vocabulary.h:275 |
| 10 | token added in 10.0 | `blind_via` | export to KiCad 9.x: not expressible (lossy) | include/drc_rules_release_vocabulary.h:276 |
| 11 | token added in 10.0 | `bridged_mask` | export to KiCad 9.x: not expressible (lossy) | include/drc_rules_release_vocabulary.h:278 |
| 12 | token added in 10.0 | `solder_mask_expansion` | export to KiCad 9.x: not expressible (lossy) | include/drc_rules_release_vocabulary.h:320 |
| 13 | token added in 10.0 | `solder_mask_sliver` | export to KiCad 9.x: not expressible (lossy) | include/drc_rules_release_vocabulary.h:321 |
| 14 | token added in 10.0 | `solder_paste_abs_margin` | export to KiCad 9.x: not expressible (lossy) | include/drc_rules_release_vocabulary.h:322 |
| 15 | token added in 10.0 | `solder_paste_rel_margin` | export to KiCad 9.x: not expressible (lossy) | include/drc_rules_release_vocabulary.h:323 |
| 16 | token added in 10.0 | `through_via` | export to KiCad 9.x: not expressible (lossy) | include/drc_rules_release_vocabulary.h:331 |
| 17 | token added in 10.0 | `via_dangling` | export to KiCad 9.x: not expressible (lossy) | include/drc_rules_release_vocabulary.h:339 |
| 18 | function added after 10.0 | `customproperty()` | export to KiCad 10.0 or older: not expressible (lossy); import: no effect | pcbnew/pcbexpr_functions.cpp:1964 |
| 19 | function added after 10.0 | `hascustomproperty()` | export to KiCad 10.0 or older: not expressible (lossy); import: no effect | pcbnew/pcbexpr_functions.cpp:1965 |
| 20 | function added after 10.0 | `hasnetchain()` | export to KiCad 10.0 or older: not expressible (lossy); import: no effect | pcbnew/pcbexpr_functions.cpp:1955 |
| 21 | function added after 10.0 | `innetchain()` | export to KiCad 10.0 or older: not expressible (lossy); import: no effect | pcbnew/pcbexpr_functions.cpp:1954 |
| 22 | function added after 10.0 | `innetchainclass()` | export to KiCad 10.0 or older: not expressible (lossy); import: no effect | pcbnew/pcbexpr_functions.cpp:1956 |
| 23 | function added after 10.0 | `intersectskeepout()` | export to KiCad 10.0 or older: not expressible (lossy); import: no effect | pcbnew/pcbexpr_functions.cpp:1936 |
| 24 | function added after 10.0 | `isstackedvia()` | export to KiCad 10.0 or older: not expressible (lossy); import: no effect | pcbnew/pcbexpr_functions.cpp:1943 |
| 25 | function added in 10.0 | `hasexactnetclass()` | export to KiCad 9.x: not expressible (lossy) | include/drc_rules_release_vocabulary.h:353 |
| 26 | function added in 10.0 | `isblindvia()` | export to KiCad 9.x: not expressible (lossy) | include/drc_rules_release_vocabulary.h:365 |
| 27 | function added in 10.0 | `isburiedvia()` | export to KiCad 9.x: not expressible (lossy) | include/drc_rules_release_vocabulary.h:366 |
| 28 | function added in 10.0 | `memberofsheetorchildren()` | export to KiCad 9.x: not expressible (lossy) | include/drc_rules_release_vocabulary.h:374 |
| 29 | property removed by 10.0 | `filled` | prop() names are version-gated; the exporter checks the target's vocabulary | include/drc_rules_release_vocabulary.h:156 |
| 30 | property removed by 10.0 | `keep_out_copper_pours` | prop() names are version-gated; the exporter checks the target's vocabulary | include/drc_rules_release_vocabulary.h:172 |
| 31 | property removed by 10.0 | `value` | prop() names are version-gated; the exporter checks the target's vocabulary | include/drc_rules_release_vocabulary.h:258 |
| 32 | property rename | `keep_out_copper_pours` (9.0) vs `keep_out_zone_fills` (10.0, master `Keep Out Zone Fills`); a rule using the old name stops resolving. | prop() names are version-gated; the exporter checks the target's vocabulary | include/drc_rules_release_vocabulary.h:172 |
| 33 | properties added in 10.0 | `auto_thickness`, `back_covering`, `back_plugging`, `backdrill_mode`, `barcode_type`, `best_length_ratio`, `best_width_ratio`, `bottom_backdrill_size`, `bottom_counterbore_depth`, `bottom_countersink_angle`, `capping` ... | recorded for the exporter (no mapping change) | include/drc_rules_release_vocabulary.h:377-553 |
| 34 | properties added after 10.0 | 83 expression names registered at this commit but absent from the 10.0 array: `affects_cursor_snap`, `affects_placement`, `affects_routing`, `angle`, `arc_end_angle`, `arc_start_angle`, `arrow_length` ... | recorded for the exporter (no mapping change) | include/drc_rules_release_vocabulary.h:377-553 |
| 35 | units | Time (ps, fs) and angle (deg) suffixes are not readable by every 9.0.x release. | export to 9.x: no deg, ps or fs units | include/drc_rules_release_vocabulary.h:265-269 |
| 36 | text variables | `${...}` in a rules file is only expanded from 10.0; the downgrade filter comments out such rules for 9.0. `${Class:...}` selectors are read as a reference designator before 10.0.6. | export to 9.x: expand variables before writing | include/drc_rules_release_vocabulary.h:575 |
| 37 | semantic break | `intersectsArea()` on a zone item tested the zone fill rather than the outline; the downgrade filter treats it as unsupported for older targets. | recorded for the exporter (no mapping change) | include/drc_rules_downgrade.h:53-57 |
| 38 | parent navigation | `A.Parent.<Field>` / `A.Parent` as an object is new; the standalone `Parent` string property is older. | recorded for the exporter (no mapping change) | include/drc_rules_downgrade.h:409-411 |
| 39 | property aliases | `Origin_X`/`Origin_Y` are aliased to `Start X`/`Start Y`; neither alias appears in the 9.0 or 10.0 property arrays. | prop() names are version-gated; the exporter checks the target's vocabulary | pcbnew/pcbexpr_evaluator.cpp:807-811 |
| 40 | downgrade export | Saving a project for an older KiCad filters the .kicad_dru: rules the target cannot parse are commented out and reported; surviving rules stay byte-identical; the version header is rewritten. | recorded for the exporter (no mapping change) | include/drc_rules_downgrade.h:35-37 |
| 41 | deprecated keywords | `mechanical_clearance` -> `physical_clearance`, `mechanical_hole_clearance` -> `physical_hole_clearance`, `hole` -> `hole_size` (parser warns). | import accepts; export writes the current spelling | pcbnew/drc/drc_rule_parser.cpp:520-534 |
| 42 | deprecated functions | `insideCourtyard`, `insideFrontCourtyard`, `insideBackCourtyard`, `insideArea`, `memberOf` are registered with a DEPRECATED tag; `insideArea` maps to the keepout implementation. | import accepts; export writes the current spelling | pcbnew/pcbexpr_functions.cpp:1927-1945 |
| 43 | dead keywords | `board_edge`, `npth`, `pth` are lexer keywords that the DRC rule parser never consumes. | recorded for the exporter (no mapping change) | common/drc_rules.keywords:4 |
| 44 | file version constant | DRC_RULE_FILE_VERSION is 20200610; the parser does not branch on the version number other than setting a too-recent flag. | recorded for the exporter (no mapping change) | pcbnew/drc/drc_rule_parser.h:34 |
| 45 | help vs parser | Help lists `(time_domain yes)`; the parser has no such token (time units select the time domain). Help omits `solder_mask_sliver` and the functions hasExactNetclass, intersectsKeepout, isBlindBuriedVia, customProperty ... | recorded for the exporter (no mapping change) | pcbnew/dialogs/panel_setup_rules_help_6notes.md:34 |
| 46 | help vs parser | The 'Missing constraint type' message omits creepage, text_height, text_thickness, connection_width, stub_length, return_path, bridged_mask, solder_mask_sliver; the switch still accepts them. | recorded for the exporter (no mapping change) | pcbnew/drc/drc_rule_parser.cpp:535-548 |
| 47 | graphical rule editor | Editor panels compile to ordinary DRU text: Via Style -> via_diameter + hole_size; permitted layers and allowed orientation -> assertion; vias under SMD -> disallow via; unmapped constraints go to a custom-rule panel. | recorded for the exporter (no mapping change) | pcbnew/dialogs/panel_setup_rules_help_6notes.md:13-25 |
| 48 | IPC API | Custom rules are also exchanged through the API as CustomRule protos; DRC_RULE::FormatRuleFromProto renders them to DRU text (DIFF_PAIR_INTRA_SKEW becomes skew + within_diff_pairs). | recorded for the exporter (no mapping change) | pcbnew/drc/drc_rule.cpp:121-197 |

#### 4.1.11 DRC test providers

| # | Provider | Canonical kinds checked | Provenance |
|---|---|---|---|
| 1 | `annular_width` | annular_width | pcbnew/drc/drc_test_provider_annular_width.cpp:31-36, :52, :464 |
| 2 | `copper width` | connection_width | pcbnew/drc/drc_test_provider_connection_width.cpp:49-54, :108, :321-324, :586 |
| 3 | `connectivity` | unrouted_net, dangling_via, dangling_track, isolated_copper, post_machining_conflict, track_not_centered_on_via | pcbnew/drc/drc_test_provider_connectivity.cpp:32-38, :50, :94-98, :332 |
| 4 | `clearance` | clearance, hole_clearance, tracks_crossing, short_circuit | pcbnew/drc/drc_test_provider_copper_clearance.cpp:48-57, :71, :128, :1513 |
| 5 | `courtyard_clearance` | courtyard_clearance, courtyard_required, hole_in_courtyard | pcbnew/drc/drc_test_provider_courtyard_clearance.cpp:28-36, :52, :300-313, :400 |
| 6 | `creepage` | creepage | pcbnew/drc/drc_test_provider_creepage.cpp:40-43, :58, :293, :307-308, :384, :441 |
| 7 | `diff_pair_coupling` | diff_pair_gap, diff_pair_uncoupled | pcbnew/drc/drc_test_provider_diff_pair_coupling.cpp:36-42, :57, :171, :777 |
| 8 | `disallow` | disallow, outline_integrity[text_on_edge_cuts] | pcbnew/drc/drc_test_provider_disallow.cpp:36-40, :53, :398 |
| 9 | `edge_clearance` | edge_clearance, silk_to_edge_clearance, outline_integrity[copper_outside_board] | pcbnew/drc/drc_test_provider_edge_clearance.cpp:41-48, :72, :176, :358, :652 |
| 10 | `footprint checks` | footprint_integrity, short_circuit | pcbnew/drc/drc_test_provider_footprint_checks.cpp:28-39, :53, :152 |
| 11 | `hole_size` | hole_size | pcbnew/drc/drc_test_provider_hole_size.cpp:30-36, :48, :261 |
| 12 | `hole_to_hole_clearance` | hole_to_hole, hole_colocated | pcbnew/drc/drc_test_provider_hole_to_hole.cpp:32-36, :52, :316 |
| 13 | `library_parity` | library_parity | pcbnew/drc/drc_test_provider_library_parity.cpp:40-44, :60, :1302 |
| 14 | `length` | length, net_chain_length, stub_length, return_path, skew, via_count | pcbnew/drc/drc_test_provider_matched_length.cpp:42-48, :62, :1031 |
| 15 | `miscellaneous` | outline_integrity[invalid_outline], item_on_disabled_layer, text_variable_unresolved, custom_assertion, tuning_profile_consistency[missing_profile] | pcbnew/drc/drc_test_provider_misc.cpp:46-56, :129, :621, :717 |
| 16 | `net chain tuning profiles` | tuning_profile_consistency[chain_mixed_profiles] | pcbnew/drc/drc_test_provider_net_chain_tuning_profiles.cpp:31-34, :45, :101 |
| 17 | `physical_clearance` | physical_clearance, physical_hole_clearance | pcbnew/drc/drc_test_provider_physical_clearance.cpp:36-41, :55, :78, :873 |
| 18 | `schematic_parity` | board_schematic_parity | pcbnew/drc/drc_test_provider_schematic_parity.cpp:31-40, :58, :409 |
| 19 | `silk_clearance` | silk_clearance | pcbnew/drc/drc_test_provider_silk_clearance.cpp:36-40, :55, :74, :333 |
| 20 | `sliver checker` | copper_sliver | pcbnew/drc/drc_test_provider_sliver_checker.cpp:34-37, :67, :88, :91, :193, :286 |
| 21 | `solder_mask_issues` | silk_to_mask_clearance, solder_mask_bridge, solder_mask_sliver | pcbnew/drc/drc_test_provider_solder_mask.cpp:43-47, :84, :517, :591, :1003-1004, :1048, :1089 |
| 22 | `text_dimensions` | text_height, text_thickness | pcbnew/drc/drc_test_provider_text_dims.cpp:33-37, :50, :285 |
| 23 | `text_mirroring` | text_mirroring | pcbnew/drc/drc_test_provider_text_mirroring.cpp:33-38, :50, :132 |
| 24 | `track angle` | track_angle | pcbnew/drc/drc_test_provider_track_angle.cpp:31-34, :47, :194 |
| 25 | `segment_length` | track_segment_length | pcbnew/drc/drc_test_provider_track_segment_length.cpp:28-31, :44, :170 |
| 26 | `width` | track_width | pcbnew/drc/drc_test_provider_track_width.cpp:27-30, :43, :163 |
| 27 | `diameter` | via_diameter | pcbnew/drc/drc_test_provider_via_diameter.cpp:26-29, :43, :191 |
| 28 | `microvia_stack` | microvia_stack, hole_aspect_ratio | pcbnew/drc/drc_test_provider_via_stack.cpp:40-68, :79, :164, :184, :251, :466 |
| 29 | `via_stitch_overlap` | via_stitch_overlap | pcbnew/drc/drc_test_provider_via_stitch.cpp:28-33, :44, :106 |
| 30 | `zone connections` | starved_thermal, zone_connection, thermal_relief_gap | pcbnew/drc/drc_test_provider_zone_connections.cpp:36-40, :54, :110-122, :370 |
| 31 | `interactive_courtyard_clearance` | courtyard_clearance (interactive use only) | pcbnew/drc/drc_interactive_courtyard_clearance.h:27-44 |

#### 4.1.12 DRC codes and severity keys (all 78)

| # | DRC code | Settings key | Default severity | Canonical kind[variant] | Quality | Provenance (enum; key) |
|---|---|---|---|---|---|---|
| 1 | `DRCE_UNCONNECTED_ITEMS` | `unconnected_items` | error | unrouted_net | exact | pcbnew/drc/drc_item.h:48; pcbnew/drc/drc_item.cpp:95 |
| 2 | `DRCE_SHORTING_ITEMS` | `shorting_items` | error | short_circuit | exact | pcbnew/drc/drc_item.h:49; pcbnew/drc/drc_item.cpp:99 |
| 3 | `DRCE_ALLOWED_ITEMS` | `items_not_allowed` | error | disallow | exact | pcbnew/drc/drc_item.h:50; pcbnew/drc/drc_item.cpp:103 |
| 4 | `DRCE_TEXT_ON_EDGECUTS` | `text_on_edge_cuts` | error | outline_integrity[text_on_edge_cuts] | exact | pcbnew/drc/drc_item.h:51; pcbnew/drc/drc_item.cpp:107 |
| 5 | `DRCE_CLEARANCE` | `clearance` | error | clearance (also physical_clearance) | exact | pcbnew/drc/drc_item.h:52; pcbnew/drc/drc_item.cpp:110 |
| 6 | `DRCE_CREEPAGE` | `creepage` | error | creepage | exact | pcbnew/drc/drc_item.h:53; pcbnew/drc/drc_item.cpp:114 |
| 7 | `DRCE_TRACKS_CROSSING` | `tracks_crossing` | error | tracks_crossing | exact | pcbnew/drc/drc_item.h:54; pcbnew/drc/drc_item.cpp:118 |
| 8 | `DRCE_EDGE_CLEARANCE` | `copper_edge_clearance` | error | edge_clearance | exact | pcbnew/drc/drc_item.h:55; pcbnew/drc/drc_item.cpp:122 |
| 9 | `DRCE_ZONES_INTERSECT` | none | n/a (no item) | none | none: enum value with no DRC item (dead code) | pcbnew/drc/drc_item.h:56; no DRC_ITEM |
| 10 | `DRCE_ISOLATED_COPPER` | `isolated_copper` | warning | isolated_copper | exact | pcbnew/drc/drc_item.h:57; pcbnew/drc/drc_item.cpp:130 |
| 11 | `DRCE_STARVED_THERMAL` | `starved_thermal` | error | starved_thermal | exact | pcbnew/drc/drc_item.h:58; pcbnew/drc/drc_item.cpp:134 |
| 12 | `DRCE_DANGLING_VIA` | `via_dangling` | warning | dangling_via | exact | pcbnew/drc/drc_item.h:59; pcbnew/drc/drc_item.cpp:138 |
| 13 | `DRCE_DANGLING_TRACK` | `track_dangling` | warning | dangling_track | exact | pcbnew/drc/drc_item.h:60; pcbnew/drc/drc_item.cpp:142 |
| 14 | `DRCE_DRILLED_HOLES_TOO_CLOSE` | `hole_to_hole` | warning | hole_to_hole | exact | pcbnew/drc/drc_item.h:61; pcbnew/drc/drc_item.cpp:150 |
| 15 | `DRCE_DRILLED_HOLES_COLOCATED` | `holes_co_located` | warning | hole_colocated | exact | pcbnew/drc/drc_item.h:62; pcbnew/drc/drc_item.cpp:154 |
| 16 | `DRCE_HOLE_CLEARANCE` | `hole_clearance` | error | hole_clearance (also physical_hole_clearance) | exact | pcbnew/drc/drc_item.h:63; pcbnew/drc/drc_item.cpp:146 |
| 17 | `DRCE_TRACK_WIDTH` | `track_width` | error | track_width | exact | pcbnew/drc/drc_item.h:64; pcbnew/drc/drc_item.cpp:162 |
| 18 | `DRCE_TRACK_ANGLE` | `track_angle` | error | track_angle | exact | pcbnew/drc/drc_item.h:65; pcbnew/drc/drc_item.cpp:166 |
| 19 | `DRCE_TRACK_SEGMENT_LENGTH` | `track_segment_length` | error | track_segment_length | exact | pcbnew/drc/drc_item.h:66; pcbnew/drc/drc_item.cpp:170 |
| 20 | `DRCE_ANNULAR_WIDTH` | `annular_width` | error | annular_width | exact | pcbnew/drc/drc_item.h:67; pcbnew/drc/drc_item.cpp:174 |
| 21 | `DRCE_CONNECTION_WIDTH` | `connection_width` | warning | connection_width | exact | pcbnew/drc/drc_item.h:68; pcbnew/drc/drc_item.cpp:158 |
| 22 | `DRCE_DRILL_OUT_OF_RANGE` | `drill_out_of_range` | error | hole_size | exact | pcbnew/drc/drc_item.h:69; pcbnew/drc/drc_item.cpp:178 |
| 23 | `DRCE_VIA_DIAMETER` | `via_diameter` | error | via_diameter | exact | pcbnew/drc/drc_item.h:70; pcbnew/drc/drc_item.cpp:182 |
| 24 | `DRCE_PADSTACK` | `padstack` | warning | footprint_integrity[padstack] | exact | pcbnew/drc/drc_item.h:71; pcbnew/drc/drc_item.cpp:186 |
| 25 | `DRCE_PADSTACK_INVALID` | `padstack_invalid` | error (fixed) | footprint_integrity[padstack_invalid] | exact | pcbnew/drc/drc_item.h:72; pcbnew/drc/drc_item.cpp:190 |
| 26 | `DRCE_MICROVIA_DRILL_OUT_OF_RANGE` | `microvia_drill_out_of_range` | error | hole_size (scope type(via.micro)) | exact | pcbnew/drc/drc_item.h:73; pcbnew/drc/drc_item.cpp:194 |
| 27 | `DRCE_MALFORMED_MICROVIA_STACK_SPAN` | `malformed_microvia_stack_span` | error | microvia_stack[malformed_span] | exact | pcbnew/drc/drc_item.h:74; pcbnew/drc/drc_item.cpp:198 |
| 28 | `DRCE_MICROVIA_STACK_NOT_FILLED` | `microvia_stack_not_filled` | warning | microvia_stack[not_filled] | exact | pcbnew/drc/drc_item.h:75; pcbnew/drc/drc_item.cpp:202 |
| 29 | `DRCE_MICROVIA_STACK_DEPTH` | `microvia_stack_depth` | error | microvia_stack[depth] | exact | pcbnew/drc/drc_item.h:76; pcbnew/drc/drc_item.cpp:205 |
| 30 | `DRCE_HOLE_ASPECT_RATIO` | `hole_aspect_ratio` | error | hole_aspect_ratio | exact | pcbnew/drc/drc_item.h:77; pcbnew/drc/drc_item.cpp:208 |
| 31 | `DRCE_MICROVIA_CROSSES_CORE` | `microvia_crosses_core` | warning | microvia_stack[crosses_core] | exact | pcbnew/drc/drc_item.h:78; pcbnew/drc/drc_item.cpp:211 |
| 32 | `DRCE_OVERLAPPING_FOOTPRINTS` | `courtyards_overlap` | error | courtyard_clearance | exact | pcbnew/drc/drc_item.h:79; pcbnew/drc/drc_item.cpp:214 |
| 33 | `DRCE_MISSING_COURTYARD` | `missing_courtyard` | ignore | courtyard_required[missing] | exact | pcbnew/drc/drc_item.h:80; pcbnew/drc/drc_item.cpp:218 |
| 34 | `DRCE_MALFORMED_COURTYARD` | `malformed_courtyard` | error | courtyard_required[malformed] | exact | pcbnew/drc/drc_item.h:81; pcbnew/drc/drc_item.cpp:222 |
| 35 | `DRCE_PTH_IN_COURTYARD` | `pth_inside_courtyard` | error | hole_in_courtyard[pth] | exact | pcbnew/drc/drc_item.h:83; pcbnew/drc/drc_item.cpp:226 |
| 36 | `DRCE_NPTH_IN_COURTYARD` | `npth_inside_courtyard` | error | hole_in_courtyard[npth] | exact | pcbnew/drc/drc_item.h:84; pcbnew/drc/drc_item.cpp:230 |
| 37 | `DRCE_DISABLED_LAYER_ITEM` | `item_on_disabled_layer` | error | item_on_disabled_layer | exact | pcbnew/drc/drc_item.h:85; pcbnew/drc/drc_item.cpp:234 |
| 38 | `DRCE_INVALID_OUTLINE` | `invalid_outline` | error | outline_integrity[invalid_outline] | exact | pcbnew/drc/drc_item.h:86; pcbnew/drc/drc_item.cpp:238 |
| 39 | `DRCE_MISSING_FOOTPRINT` | `missing_footprint` | warning | board_schematic_parity[missing_footprint] | exact | pcbnew/drc/drc_item.h:88; pcbnew/drc/drc_item.cpp:246 |
| 40 | `DRCE_DUPLICATE_FOOTPRINT` | `duplicate_footprints` | warning | board_schematic_parity[duplicate_footprint] | exact | pcbnew/drc/drc_item.h:89; pcbnew/drc/drc_item.cpp:242 |
| 41 | `DRCE_EXTRA_FOOTPRINT` | `extra_footprint` | warning | board_schematic_parity[extra_footprint] | exact | pcbnew/drc/drc_item.h:90; pcbnew/drc/drc_item.cpp:250 |
| 42 | `DRCE_NET_CONFLICT` | `net_conflict` | warning | board_schematic_parity[net_conflict] | exact | pcbnew/drc/drc_item.h:91; pcbnew/drc/drc_item.cpp:254 |
| 43 | `DRCE_SCHEMATIC_PARITY` | `footprint_symbol_mismatch` | warning | board_schematic_parity[symbol_mismatch] | exact | pcbnew/drc/drc_item.h:92; pcbnew/drc/drc_item.cpp:258 |
| 44 | `DRCE_FOOTPRINT_FILTERS` | `footprint_filters_mismatch` | ignore | board_schematic_parity[footprint_filter_mismatch] | exact | pcbnew/drc/drc_item.h:93; pcbnew/drc/drc_item.cpp:262 |
| 45 | `DRCE_FOOTPRINT_TYPE_MISMATCH` | `footprint_type_mismatch` | ignore | footprint_integrity[type_mismatch] | exact | pcbnew/drc/drc_item.h:95; pcbnew/drc/drc_item.cpp:363 |
| 46 | `DRCE_LIB_FOOTPRINT_ISSUES` | `lib_footprint_issues` | warning | library_parity[missing] | exact | pcbnew/drc/drc_item.h:96; pcbnew/drc/drc_item.cpp:274 |
| 47 | `DRCE_LIB_FOOTPRINT_MISMATCH` | `lib_footprint_mismatch` | warning | library_parity[mismatch] | exact | pcbnew/drc/drc_item.h:97; pcbnew/drc/drc_item.cpp:278 |
| 48 | `DRCE_PAD_TH_WITH_NO_HOLE` | `through_hole_pad_without_hole` | error | footprint_integrity[th_pad_without_hole] | exact | pcbnew/drc/drc_item.h:98; pcbnew/drc/drc_item.cpp:367 |
| 49 | `DRCE_FOOTPRINT` | `footprint` | error | footprint_integrity[net_tie_definition] | exact | pcbnew/drc/drc_item.h:99; pcbnew/drc/drc_item.cpp:359 |
| 50 | `DRCE_FOOTPRINT_SCALED_WITH_PADS` | `footprint_scaled_with_pads` | warning | footprint_integrity[scaled_with_pads] | exact | pcbnew/drc/drc_item.h:100; pcbnew/drc/drc_item.cpp:371 |
| 51 | `DRCE_UNRESOLVED_VARIABLE` | `unresolved_variable` | error | text_variable_unresolved | exact | pcbnew/drc/drc_item.h:102; pcbnew/drc/drc_item.cpp:282 |
| 52 | `DRCE_ASSERTION_FAILURE` | `assertion_failure` | error | custom_assertion | exact | pcbnew/drc/drc_item.h:103; pcbnew/drc/drc_item.cpp:286 |
| 53 | `DRCE_GENERIC_WARNING` | `generic_warning` | error in the severity map (loop default) | custom_assertion (text warnings; fixed severity) | exact | pcbnew/drc/drc_item.h:104; pcbnew/drc/drc_item.cpp:290 |
| 54 | `DRCE_GENERIC_ERROR` | `generic_error` | error (fixed) | custom_assertion (text errors; fixed severity) | exact | pcbnew/drc/drc_item.h:105; pcbnew/drc/drc_item.cpp:294 |
| 55 | `DRCE_COPPER_SLIVER` | `copper_sliver` | warning | copper_sliver | exact | pcbnew/drc/drc_item.h:107; pcbnew/drc/drc_item.cpp:298 |
| 56 | `DRCE_SOLDERMASK_BRIDGE` | `solder_mask_bridge` | error | solder_mask_bridge | exact | pcbnew/drc/drc_item.h:108; pcbnew/drc/drc_item.cpp:302 |
| 57 | `DRCE_SILK_MASK_CLEARANCE` | `silk_over_copper` | warning | silk_to_mask_clearance | exact | pcbnew/drc/drc_item.h:111; pcbnew/drc/drc_item.cpp:306 |
| 58 | `DRCE_SILK_EDGE_CLEARANCE` | `silk_edge_clearance` | warning | silk_to_edge_clearance | exact | pcbnew/drc/drc_item.h:113; pcbnew/drc/drc_item.cpp:310 |
| 59 | `DRCE_SILK_CLEARANCE` | `silk_overlap` | warning | silk_clearance | exact | pcbnew/drc/drc_item.h:114; pcbnew/drc/drc_item.cpp:314 |
| 60 | `DRCE_TEXT_HEIGHT` | `text_height` | warning | text_height | exact | pcbnew/drc/drc_item.h:115; pcbnew/drc/drc_item.cpp:318 |
| 61 | `DRCE_TEXT_THICKNESS` | `text_thickness` | warning | text_thickness | exact | pcbnew/drc/drc_item.h:116; pcbnew/drc/drc_item.cpp:322 |
| 62 | `DRCE_LENGTH_OUT_OF_RANGE` | `length_out_of_range` | error | length (also net_chain_length) | exact | pcbnew/drc/drc_item.h:118; pcbnew/drc/drc_item.cpp:326 |
| 63 | `DRCE_NET_CHAIN_STUB_TOO_LONG` | `net_chain_stub_length` | error | stub_length | exact | pcbnew/drc/drc_item.h:119; pcbnew/drc/drc_item.cpp:330 |
| 64 | `DRCE_NET_CHAIN_RETURN_PATH_BREAK` | `net_chain_return_path` | error | return_path | exact | pcbnew/drc/drc_item.h:120; pcbnew/drc/drc_item.cpp:334 |
| 65 | `DRCE_NET_CHAIN_TUNING_PROFILES` | `net_chain_tuning_profiles` | error | tuning_profile_consistency[chain_mixed_profiles] | exact | pcbnew/drc/drc_item.h:121; pcbnew/drc/drc_item.cpp:338 |
| 66 | `DRCE_SKEW_OUT_OF_RANGE` | `skew_out_of_range` | error | skew | exact | pcbnew/drc/drc_item.h:122; pcbnew/drc/drc_item.cpp:342 |
| 67 | `DRCE_VIA_COUNT_OUT_OF_RANGE` | `too_many_vias` | error | via_count | exact | pcbnew/drc/drc_item.h:123; pcbnew/drc/drc_item.cpp:347 |
| 68 | `DRCE_DP_GAP_OUT_OF_RANGE` | `diff_pair_gap_out_of_range` | error | diff_pair_gap | exact | pcbnew/drc/drc_item.h:124; pcbnew/drc/drc_item.cpp:351 |
| 69 | `DRCE_DP_UNCOUPLED_LENGTH_TOO_LONG` | `diff_pair_uncoupled_length_too_long` | error | diff_pair_uncoupled | exact | pcbnew/drc/drc_item.h:125; pcbnew/drc/drc_item.cpp:355 |
| 70 | `DRCE_MIRRORED_TEXT_ON_FRONT_LAYER` | `mirrored_text_on_front_layer` | warning | text_mirroring[front_mirrored] | exact | pcbnew/drc/drc_item.h:127; pcbnew/drc/drc_item.cpp:375 |
| 71 | `DRCE_UNMIRRORED_TEXT_ON_BACK_LAYER` | `nonmirrored_text_on_back_layer` | warning | text_mirroring[back_unmirrored] | exact | pcbnew/drc/drc_item.h:128; pcbnew/drc/drc_item.cpp:379 |
| 72 | `DRCE_MISSING_TUNING_PROFILE` | `missing_tuning_profile` | warning | tuning_profile_consistency[missing_profile] | exact | pcbnew/drc/drc_item.h:130; pcbnew/drc/drc_item.cpp:383 |
| 73 | `DRCE_TUNING_PROFILE_IMPLICIT_RULES` | `tuning_profile_track_geometries` | ignore | tuning_profile_consistency[implicit_geometry] (severity carrier, never emitted) | exact | pcbnew/drc/drc_item.h:131; pcbnew/drc/drc_item.cpp:387 |
| 74 | `DRCE_TRACK_ON_POST_MACHINED_LAYER` | `track_on_post_machined_layer` | error | post_machining_conflict | exact | pcbnew/drc/drc_item.h:133; pcbnew/drc/drc_item.cpp:391 |
| 75 | `DRCE_TRACK_NOT_CENTERED_ON_VIA` | `track_not_centered_on_via` | ignore | track_not_centered_on_via | exact | pcbnew/drc/drc_item.h:135; pcbnew/drc/drc_item.cpp:395 |
| 76 | `DRCE_SCHEMATIC_FIELDS_PARITY` | `footprint_symbol_field_mismatch` | warning | board_schematic_parity[field_mismatch] | exact | pcbnew/drc/drc_item.h:137; pcbnew/drc/drc_item.cpp:266 |
| 77 | `DRCE_VIA_STITCH_OVERLAP` | `via_stitch_overlap` | error | via_stitch_overlap | exact | pcbnew/drc/drc_item.h:139; pcbnew/drc/drc_item.cpp:270 |
| 78 | `DRCE_COPPER_OUTSIDE_BOARD` | `copper_outside_board` | warning | outline_integrity[copper_outside_board] | exact | pcbnew/drc/drc_item.h:142; pcbnew/drc/drc_item.cpp:126 |

`DRCE_FIRST` and `DRCE_LAST` are range aliases, not codes: `DRCE_FIRST` = 1 = `DRCE_UNCONNECTED_ITEMS` and
`DRCE_LAST` = `DRCE_COPPER_OUTSIDE_BOARD` (pcbnew/drc/drc_item.h:47-48, :144). The loop over that range at
pcbnew/board_design_settings.cpp:169-170 sets `error` as the default for every code, and lines :172-219 then override
the codes whose default above is not `error`.

#### 4.1.13 ERC codes and severity keys (all 63)

| # | ERC code | Settings key | Default severity | Canonical kind[variant] | Quality | Provenance |
|---|---|---|---|---|---|---|
| 1 | `ERCE_DUPLICATE_SHEET_NAME` | `duplicate_sheet_names` | error | erc_annotation[duplicate_sheet_name] | exact | eeschema/erc/erc_settings.h:43; erc_item.cpp#L42 |
| 2 | `ERCE_ENDPOINT_OFF_GRID` | `endpoint_off_grid` | warning | erc_drawing[endpoint_off_grid] | exact | eeschema/erc/erc_settings.h:44; erc_item.cpp#L46 |
| 3 | `ERCE_PIN_NOT_CONNECTED` | `pin_not_connected` | error | erc_unconnected_pin | exact | eeschema/erc/erc_settings.h:45; erc_item.cpp#L50 |
| 4 | `ERCE_PIN_NOT_DRIVEN` | `pin_not_driven` | error | erc_undriven_input[input] | exact | eeschema/erc/erc_settings.h:46; erc_item.cpp#L54 |
| 5 | `ERCE_POWERPIN_NOT_DRIVEN` | `power_pin_not_driven` | error | erc_undriven_input[power_input] | exact | eeschema/erc/erc_settings.h:48; erc_item.cpp#L58 |
| 6 | `ERCE_HIERACHICAL_LABEL` | `hier_label_mismatch` | error | erc_label[hier_label_mismatch] | exact | eeschema/erc/erc_settings.h:50; erc_item.cpp#L86 |
| 7 | `ERCE_NOCONNECT_CONNECTED` | `no_connect_connected` | warning | erc_no_connect[connected] | exact | eeschema/erc/erc_settings.h:51; erc_item.cpp#L98 |
| 8 | `ERCE_NOCONNECT_NOT_CONNECTED` | `no_connect_dangling` | warning | erc_no_connect[dangling] | exact | eeschema/erc/erc_settings.h:52; erc_item.cpp#L102 |
| 9 | `ERCE_LABEL_NOT_CONNECTED` | `label_dangling` | error | erc_dangling[label] | exact | eeschema/erc/erc_settings.h:53; erc_item.cpp#L106 |
| 10 | `ERCE_SIMILAR_LABELS` | `similar_labels` | warning | erc_label[similar_labels] | exact | eeschema/erc/erc_settings.h:54; erc_item.cpp#L114 |
| 11 | `ERCE_SIMILAR_POWER` | `similar_power` | warning | erc_label[similar_power] | exact | eeschema/erc/erc_settings.h:55; erc_item.cpp#L118 |
| 12 | `ERCE_SIMILAR_LABEL_AND_POWER` | `similar_label_and_power` | warning | erc_label[similar_label_and_power] | exact | eeschema/erc/erc_settings.h:56; erc_item.cpp#L122 |
| 13 | `ERCE_SINGLE_GLOBAL_LABEL` | `single_global_label` | ignore | erc_label[single_global_label] | exact | eeschema/erc/erc_settings.h:57; erc_item.cpp#L126 |
| 14 | `ERCE_SAME_LOCAL_GLOBAL_LABEL` | `same_local_global_label` | warning | erc_label[same_local_global_label] | exact | eeschema/erc/erc_settings.h:58; erc_item.cpp#L130 |
| 15 | `ERCE_SAME_LOCAL_GLOBAL_POWER` | `same_local_global_power` | warning | erc_label[same_local_global_power] | exact | eeschema/erc/erc_settings.h:59; erc_item.cpp#L134 |
| 16 | `ERCE_DIFFERENT_UNIT_FP` | `different_unit_footprint` | error | erc_multi_unit[different_unit_footprint] | exact | eeschema/erc/erc_settings.h:60; erc_item.cpp#L138 |
| 17 | `ERCE_MISSING_POWER_INPUT_PIN` | `missing_power_pin` | error | erc_multi_unit[missing_power_pin] | exact | eeschema/erc/erc_settings.h:62; erc_item.cpp#L260 |
| 18 | `ERCE_MISSING_INPUT_PIN` | `missing_input_pin` | warning | erc_multi_unit[missing_input_pin] | exact | eeschema/erc/erc_settings.h:64; erc_item.cpp#L252 |
| 19 | `ERCE_MISSING_BIDI_PIN` | `missing_bidi_pin` | warning | erc_multi_unit[missing_bidi_pin] | exact | eeschema/erc/erc_settings.h:65; erc_item.cpp#L256 |
| 20 | `ERCE_MISSING_UNIT` | `missing_unit` | warning | erc_multi_unit[missing_unit] | exact | eeschema/erc/erc_settings.h:66; erc_item.cpp#L248 |
| 21 | `ERCE_DIFFERENT_UNIT_NET` | `different_unit_net` | error | erc_multi_unit[different_unit_net] | exact | eeschema/erc/erc_settings.h:67; erc_item.cpp#L142 |
| 22 | `ERCE_BUS_ALIAS_CONFLICT` | `bus_definition_conflict` | error | erc_bus[definition_conflict] | exact | eeschema/erc/erc_settings.h:69; erc_item.cpp#L146 |
| 23 | `ERCE_DRIVER_CONFLICT` | `multiple_net_names` | warning | erc_label[multiple_net_names] | exact | eeschema/erc/erc_settings.h:70; erc_item.cpp#L150 |
| 24 | `ERCE_BUS_ENTRY_CONFLICT` | `net_not_bus_member` | warning | erc_bus[net_not_bus_member] | exact | eeschema/erc/erc_settings.h:71; erc_item.cpp#L154 |
| 25 | `ERCE_BUS_TO_BUS_CONFLICT` | `bus_to_bus_conflict` | error | erc_bus[bus_to_bus] | exact | eeschema/erc/erc_settings.h:72; erc_item.cpp#L158 |
| 26 | `ERCE_BUS_TO_NET_CONFLICT` | `bus_to_net_conflict` | error | erc_bus[bus_to_net] | exact | eeschema/erc/erc_settings.h:74; erc_item.cpp#L162 |
| 27 | `ERCE_GROUND_PIN_NOT_GROUND` | `ground_pin_not_ground` | warning | erc_power[ground_pin_not_ground] | exact | eeschema/erc/erc_settings.h:76; erc_item.cpp#L166 |
| 28 | `ERCE_LABEL_SINGLE_PIN` | `isolated_pin_label` | warning | erc_label[isolated_pin_label] | exact | eeschema/erc/erc_settings.h:77; erc_item.cpp#L110 |
| 29 | `ERCE_UNRESOLVED_VARIABLE` | `unresolved_variable` | error | text_variable_unresolved | exact | eeschema/erc/erc_settings.h:78; erc_item.cpp#L208 |
| 30 | `ERCE_UNDEFINED_NETCLASS` | `undefined_netclass` | error | erc_netclass | exact | eeschema/erc/erc_settings.h:79; erc_item.cpp#L212 |
| 31 | `ERCE_SIMULATION_MODEL` | `simulation_model_issue` | ignore | erc_simulation_model | exact | eeschema/erc/erc_settings.h:80; erc_item.cpp#L216 |
| 32 | `ERCE_WIRE_DANGLING` | `wire_dangling` | error | erc_dangling[wire] | exact | eeschema/erc/erc_settings.h:81; erc_item.cpp#L220 |
| 33 | `ERCE_LIB_SYMBOL_ISSUES` | `lib_symbol_issues` | warning | erc_symbol_library[lib_symbol_issues] | exact | eeschema/erc/erc_settings.h:82; erc_item.cpp#L224 |
| 34 | `ERCE_LIB_SYMBOL_MISMATCH` | `lib_symbol_mismatch` | warning | erc_symbol_library[lib_symbol_mismatch] | exact | eeschema/erc/erc_settings.h:83; erc_item.cpp#L228 |
| 35 | `ERCE_FOOTPRINT_LINK_ISSUES` | `footprint_link_issues` | warning | erc_symbol_library[footprint_link_issues] | exact | eeschema/erc/erc_settings.h:84; erc_item.cpp#L232 |
| 36 | `ERCE_FOOTPRINT_FILTERS` | `footprint_filter` | ignore | erc_symbol_library[footprint_filter] | exact | eeschema/erc/erc_settings.h:86; erc_item.cpp#L236 |
| 37 | `ERCE_UNANNOTATED` | `unannotated` | error | erc_annotation[unannotated] | exact | eeschema/erc/erc_settings.h:87; erc_item.cpp#L240 |
| 38 | `ERCE_EXTRA_UNITS` | `extra_units` | error | erc_annotation[extra_units] | exact | eeschema/erc/erc_settings.h:88; erc_item.cpp#L244 |
| 39 | `ERCE_DIFFERENT_UNIT_VALUE` | `unit_value_mismatch` | error | erc_annotation[unit_value_mismatch] | exact | eeschema/erc/erc_settings.h:89; erc_item.cpp#L264 |
| 40 | `ERCE_DUPLICATE_REFERENCE` | `duplicate_reference` | error | erc_annotation[duplicate_reference] | exact | eeschema/erc/erc_settings.h:90; erc_item.cpp#L268 |
| 41 | `ERCE_BUS_ENTRY_NEEDED` | `bus_entry_needed` | error | erc_bus[entry_needed] | exact | eeschema/erc/erc_settings.h:91; erc_item.cpp#L272 |
| 42 | `ERCE_FOUR_WAY_JUNCTION` | `four_way_junction` | ignore | erc_drawing[four_way_junction] | exact | eeschema/erc/erc_settings.h:92; erc_item.cpp#L90 |
| 43 | `ERCE_LABEL_MULTIPLE_WIRES` | `label_multiple_wires` | warning | erc_label[label_multiple_wires] | exact | eeschema/erc/erc_settings.h:93; erc_item.cpp#L94 |
| 44 | `ERCE_UNCONNECTED_WIRE_ENDPOINT` | `unconnected_wire_endpoint` | warning | erc_dangling[wire_endpoint] | exact | eeschema/erc/erc_settings.h:94; erc_item.cpp#L276 |
| 45 | `ERCE_STACKED_PIN_SYNTAX` | `stacked_pin_name` | warning | erc_symbol_library[stacked_pin_name] | exact | eeschema/erc/erc_settings.h:95; erc_item.cpp#L170 |
| 46 | `ERCE_FIELD_NAME_WHITESPACE` | `field_name_whitespace` | warning | erc_symbol_library[field_name_whitespace] | exact | eeschema/erc/erc_settings.h:96; erc_item.cpp#L174 |
| 47 | `ERCE_PIN_MAP_BAD_PAD` | `pin_map_bad_pad` | error | erc_symbol_library[pin_map_bad_pad] | exact | eeschema/erc/erc_settings.h:97; erc_item.cpp#L192 |
| 48 | `ERCE_PIN_MAP_UNMAPPED_PIN` | `pin_map_unmapped_pin` | warning | erc_symbol_library[pin_map_unmapped_pin] | exact | eeschema/erc/erc_settings.h:98; erc_item.cpp#L196 |
| 49 | `ERCE_PIN_MAP_DUPLICATE_PAD` | `pin_map_duplicate_pad` | error | erc_symbol_library[pin_map_duplicate_pad] | exact | eeschema/erc/erc_settings.h:99; erc_item.cpp#L200 |
| 50 | `ERCE_PIN_MAP_STALE_PIN` | `pin_map_stale_pin` | warning | erc_symbol_library[pin_map_stale_pin] | exact | eeschema/erc/erc_settings.h:100; erc_item.cpp#L204 |
| 51 | `ERCE_EMPTY_LABEL_NAME` | `empty_label_name` | error | erc_label[empty_label_name] | exact | eeschema/erc/erc_settings.h:101; erc_item.cpp#L178 |
| 52 | `ERCE_VARIANT_SYMBOL_INVALID` | `variant_symbol_invalid` | error | erc_symbol_library[variant_symbol_invalid] | exact | eeschema/erc/erc_settings.h:102; erc_item.cpp#L280 |
| 53 | `ERCE_VARIANT_SYMBOL_INCOMPATIBLE` | `variant_symbol_incompatible` | warning | erc_symbol_library[variant_symbol_incompatible] | exact | eeschema/erc/erc_settings.h:103; erc_item.cpp#L284 |
| 54 | `ERCE_WIRED_IMPLICIT_POWER` | `wired_implicit_power` | warning | erc_power[wired_implicit_power] | exact | eeschema/erc/erc_settings.h:104; erc_item.cpp#L66 |
| 55 | `ERCE_EMPTY_POWER_VALUE` | `empty_power_value` | error | erc_power[empty_power_value] | exact | eeschema/erc/erc_settings.h:105; erc_item.cpp#L181 |
| 56 | `ERCE_DUPLICATE_SHEET_SUFFIX` | `duplicate_sheet_suffix` | error | erc_annotation[duplicate_sheet_suffix] | exact | eeschema/erc/erc_settings.h:106; erc_item.cpp#L184 |
| 57 | `ERCE_REF_SUFFIX_MISMATCH` | `ref_suffix_mismatch` | warning | erc_annotation[ref_suffix_mismatch] | exact | eeschema/erc/erc_settings.h:107; erc_item.cpp#L188 |
| 58 | `ERCE_DUPLICATE_PIN_ERROR` | `duplicate_pins` | error (forced) | erc_multi_unit[duplicate_pins] | exact | eeschema/erc/erc_settings.h:111; erc_item.cpp#L62 |
| 59 | `ERCE_PIN_TO_PIN_WARNING` | `pin_to_pin` | warning (fixed; only ignore is | erc_pin_type_conflict (matrix cell = warning) | exact | eeschema/erc/erc_settings.h:112; erc_item.cpp#L70 |
| 60 | `ERCE_PIN_TO_PIN_ERROR` | `pin_to_pin` | error (fixed; only ignore is c | erc_pin_type_conflict (matrix cell = error) | exact | eeschema/erc/erc_settings.h:113; erc_item.cpp#L74 |
| 61 | `ERCE_ANNOTATION_ACTION` | none | n/a | none | none: an annotation action passed through the error handler, not a check | eeschema/erc/erc_settings.h:114 |
| 62 | `ERCE_GENERIC_WARNING` | `generic-warning` | warning (fixed) | custom_assertion (schematic text assertion; fixed severity) | exact | eeschema/erc/erc_settings.h:116; erc_item.cpp#L78 |
| 63 | `ERCE_GENERIC_ERROR` | `generic-error` | error (fixed) | custom_assertion (schematic text assertion; fixed severity) | exact | eeschema/erc/erc_settings.h:117; erc_item.cpp#L82 |

#### 4.1.14 ERC pin-type matrix and drive rules

The 12 x 12 pin matrix (rows I, O, Bi, 3S, Pas, NIC, UnS, PwrI, PwrO, OC, OE, NC; eeschema/erc/erc_settings.cpp:45-57) maps cell by cell to
`erc_pin_type_conflict.matrix` (values ok, warning, error). Exact. The canonical matrix uses the union of KiCad pin
types and Altium port and sheet-entry types, so Altium's 17 x 17 matrix fits once its labels are confirmed (Q-A6).

| # | Item | Canonical target | Quality | Note | Provenance |
|---|---|---|---|---|---|
| 1 | Pins that need a driver | erc_undriven_input (pins needing a driver: input, power input) | exact |  | eeschema/erc/erc.cpp:327-332 |
| 2 | Drivers for ordinary inputs | erc_undriven_input[input] driver set | exact |  | eeschema/erc/erc.cpp:307-317 |
| 3 | Drivers for power inputs | erc_undriven_input[power_input] driver set | exact |  | eeschema/erc/erc.cpp:319-325 |
| 4 | Undriven net | erc_undriven_input variant selection | exact |  | eeschema/erc/erc.cpp:1797, :2114-2116 |
| 5 | Legacy minimum-drive table | erc_undriven_input (legacy drive table) | exact | kept as native data | eeschema/erc/erc_settings.h:136-141; eeschema/erc/erc_settings.cpp:61-86 |
| 6 | `erc.pin_map` in the project | erc_pin_type_conflict.matrix (0 ok, 1 warning, 2 error, 3 unconnected) | exact |  | eeschema/erc/erc_settings.h:120-127; eeschema/erc/erc_settings.cpp:203-247 |
| 7 | Severity of matrix hits | @severity of erc_pin_type_conflict (ignore only) | model |  | eeschema/erc/erc_settings.cpp:322-341 |
| 8 | Marker placement for multi-pin conflicts | none | none: marker placement heuristic (UI) |  | eeschema/erc/erc_settings.cpp:251-265 |

#### 4.1.15 KiCad 4 and 5 board-file `(setup)` keys (golden boards with file versions 4 and 20160815-20171130)

File version 4 occurs once, in speeduino-v04 `v0.4/Interface boards/Screw terminals/Screw terminal interface.kicad_pcb`
(golden-rules/speeduino-v04.txt:14-15); it is the only golden board with `hpglpenoverlay`. Rows 1-28 cover the design
keys, rows 29-33 the plot and display keys stored in the same `(setup)` block; together they cover all 58 distinct
`(setup)` tokens of the golden KiCad 4 and 5 boards.

| # | (setup) key | Example | Unit | Canonical target | Quality | Note | Provenance |
|---|---|---|---|---|---|---|---|
| 1 | `trace_clearance` | 0.2 | mm | clearance.min [netclass Default] | lossy: mapping unverified (GK-K) |  | golden-rules/eez-psu-h24005.txt:13; equivalent: unverified (GK-K) |
| 2 | `zone_clearance` | 0.6 | mm | clearance [item_override default for zones] | lossy: mapping unverified (GK-K) |  | golden-rules/eez-psu-h24005.txt:14; equivalent: unverified (GK-K) |
| 3 | `zone_45_only` | (key only) | bool | none | none: zone drawing option |  | golden-rules/eez-psu-h24005.txt:3; equivalent: unverified (GK-K) |
| 4 | `trace_min` | 0.2 | mm | track_width.min [board_floor] | lossy: mapping unverified (GK-K) |  | golden-rules/eez-psu-h24005.txt:15; equivalent: unverified (GK-K) |
| 5 | `last_trace_width` | 0.25 | mm | none | none: router state |  | golden-rules/eez-psu-h24005.txt:4; equivalent: unverified (GK-K) |
| 6 | `user_trace_width (repeated)` | 0.3 | mm | router_presets.widths | lossy: mapping unverified (GK-K) |  | golden-rules/eez-psu-h24005.txt:5; equivalent: unverified (GK-K) |
| 7 | `segment_width` | 1 | mm | none | none: graphic line default |  | golden-rules/eez-psu-h24005.txt:16; equivalent: unverified (GK-K) |
| 8 | `edge_width` | 0.15 | mm | none | none: outline line default |  | golden-rules/eez-psu-h24005.txt:17; equivalent: unverified (GK-K) |
| 9 | `via_size` | 1 | mm | via_diameter.opt [netclass Default] | lossy: mapping unverified (GK-K) |  | golden-rules/eez-psu-h24005.txt:18; equivalent: unverified (GK-K) |
| 10 | `via_drill` | 0.6 | mm | hole_size.opt [netclass Default, vias] | lossy: mapping unverified (GK-K) |  | golden-rules/eez-psu-h24005.txt:19; equivalent: unverified (GK-K) |
| 11 | `via_min_size` | 0.4 | mm | via_diameter.min [board_floor] | lossy: mapping unverified (GK-K) |  | golden-rules/eez-psu-h24005.txt:20; equivalent: unverified (GK-K) |
| 12 | `via_min_drill` | 0.3 | mm | hole_size.min [board_floor] | lossy: mapping unverified (GK-K) |  | golden-rules/eez-psu-h24005.txt:21; equivalent: unverified (GK-K) |
| 13 | `user_via (repeated)` | 0.45 0.3 | mm diameter, mm drill | router_presets.vias | lossy: mapping unverified (GK-K) |  | golden-rules/esp32-evb.txt:123; equivalent: unverified (GK-K) |
| 14 | `uvia_size` | 0.3 | mm | via_diameter.opt [netclass Default, micro] | lossy: mapping unverified (GK-K) |  | golden-rules/eez-psu-h24005.txt:22; equivalent: unverified (GK-K) |
| 15 | `uvia_drill` | 0.1 | mm | hole_size.opt [netclass Default, micro] | lossy: mapping unverified (GK-K) |  | golden-rules/eez-psu-h24005.txt:23; equivalent: unverified (GK-K) |
| 16 | `uvias_allowed` | no | bool | via_span.allowed_types (micro) | lossy: mapping unverified (GK-K) |  | golden-rules/eez-psu-h24005.txt:24; equivalent: unverified (GK-K) |
| 17 | `uvia_min_size` | 0.2 | mm | via_diameter.min [board_floor, micro] | lossy: mapping unverified (GK-K) |  | golden-rules/eez-psu-h24005.txt:25; equivalent: unverified (GK-K) |
| 18 | `uvia_min_drill` | 0.1 | mm | hole_size.min [board_floor, micro] | lossy: mapping unverified (GK-K) |  | golden-rules/eez-psu-h24005.txt:26; equivalent: unverified (GK-K) |
| 19 | `pcb_text_width` | 0.3 | mm | none | none: copper text creation default |  | golden-rules/eez-psu-h24005.txt:27; equivalent: unverified (GK-K) |
| 20 | `pcb_text_size` | (key only) | mm | none | none: copper text creation default |  | golden-rules/eez-psu-h24005.txt:3; equivalent: unverified (GK-K) |
| 21 | `mod_edge_width` | 0.15 | mm | none | none: footprint line creation default |  | golden-rules/eez-psu-h24005.txt:28; equivalent: unverified (GK-K) |
| 22 | `mod_text_width` | 0.15 | mm | none | none: footprint text creation default |  | golden-rules/eez-psu-h24005.txt:29; equivalent: unverified (GK-K) |
| 23 | `mod_text_size` | (key only) | mm | none | none: footprint text creation default |  | golden-rules/eez-psu-h24005.txt:3; equivalent: unverified (GK-K) |
| 24 | `pad_size` | (key only) | mm | none | none: new-pad default |  | golden-rules/eez-psu-h24005.txt:3; equivalent: unverified (GK-K) |
| 25 | `pad_drill` | 0 | mm | none | none: new-pad default |  | golden-rules/eez-psu-h24005.txt:30; equivalent: unverified (GK-K) |
| 26 | `pad_to_mask_clearance` | 0.05 | mm | solder_mask_expansion.opt [board default] | lossy: mapping unverified (GK-K) |  | golden-rules/eez-psu-h24005.txt:31; equivalent: unverified (GK-K) |
| 27 | `solder_mask_min_width` | 0.09 | mm | solder_mask_bridge.min_web + solder_mask_sliver.min [board_floor], as section 4.1.4 row 2 | lossy: KiCad enforces it only as the bridge criterion (section 4.1.4 row 2) | read into `m_SolderMaskMinWidth` by the current parser (pcbnew/pcb_io/kicad_sexpr/pcb_io_kicad_sexpr_parser.cpp:3302-3305†) | golden-rules/glasgow.txt:135; equivalent: unverified (GK-K) |
| 28 | `aux_axis_origin, grid_origin` | (key only) | mm coordinates | none | none: coordinate origins |  | golden-rules/eez-psu-h24005.txt:3; golden-rules/glasgow.txt:43 |
| 29 | `pcbplotparams` and its children `creategerberjobfile`, `excludeedgelayer`, `hpglpendiameter`, `hpglpennumber`, `hpglpenoverlay`, `hpglpenspeed`, `layerselection`, `linewidth`, `mirror`, `mode`, `outputdirectory`, `outputformat`, `padsonsilk`, `plotframeref`, `plotinvisibletext`, `plotreference`, `plotvalue`, `psa4output`, `psnegative`, `scaleselection`, `useauxorigin`, `usegerberadvancedattributes`, `usegerberattributes`, `usegerberextensions` (24 children) | (block) | - | none | none: plot setting, not a rule | the current loader reads them into the plot options or skips them: HPGL `mode`, `hpglpennumber`, `hpglpenspeed`, `hpglpenoverlay` are read and discarded; `linewidth`, `padsonsilk`, `plotreference`, `plotvalue` are not keywords any more and are skipped | golden-rules/eez-psu-h24005.txt:3, :32-33; golden-rules/speeduino-v04.txt:15 (`hpglpenoverlay`); pcbnew/pcb_io/kicad_sexpr/pcb_io_kicad_sexpr_parser.cpp:3410-3431†; pcbnew/pcb_plot_params.cpp:551-828† (HPGL :719-725, default skip :819-822); common/pcb_plot_params.keywords:1-55† |
| 30 | `viasonmask` (child of `pcbplotparams`) | false | bool | via_protection.tent_top = via_protection.tent_bottom = not viasonmask [board default] | exact | a legacy plot option that the current loader converts to board via tenting | golden-rules/eez-psu-h24005.txt:34; pcbnew/pcb_plot_params.cpp:711-713†; pcbnew/pcb_io/kicad_sexpr/pcb_io_kicad_sexpr_parser.cpp:3423-3428† |
| 31 | `subtractmaskfromsilk` (child of `pcbplotparams`) | false | bool | none (value kept in `params` as a fabrication-output hint) | none: plot setting, but it changes the fabricated silkscreen (silk clipped at mask openings), so the exporter writes it back | | golden-rules/eez-psu-h24005.txt:35; pcbnew/pcb_plot_params.cpp:792-794† |
| 32 | `drillshape` (child of `pcbplotparams`) | 1 | enum 0-2 | none | none: drill-mark style of plots (drill map), not a rule | | golden-rules/eez-psu-h24005.txt:36; pcbnew/pcb_plot_params.cpp:806-808† |
| 33 | `visible_elements` | (key only) | hex bitmask | none | none: display visibility, stored in the board before KiCad 6 | | golden-rules/eez-psu-h24005.txt:3; pcbnew/pcb_io/kicad_sexpr/pcb_io_kicad_sexpr_parser.cpp:3383-3396† |

The current s-expression loader (fetched for rows 27 and 29-33) also reads the design keys of rows 1-26 in its
`parseSetup` (pcbnew/pcb_io/kicad_sexpr/pcb_io_kicad_sexpr_parser.cpp:3033-3500†); rows 1-26 keep their GK-K status until they are re-checked against it (Q-K12).
One finding from that read: `uvias_allowed` and `blind_buried_vias_allowed` are parsed and their values discarded
(pcbnew/pcb_io/kicad_sexpr/pcb_io_kicad_sexpr_parser.cpp:3190-3200†), so current KiCad ignores row 16.

#### 4.1.16 KiCad 5 `(net_class ...)` keys

| # | (net_class) key | Canonical target | Quality | Provenance |
|---|---|---|---|---|
| 1 | `clearance` | clearance.min [netclass] | lossy: legacy board parser not in checkout (GK-K) | golden-rules/eez-psu-h24005.txt:37-38; golden-rules/esp32-evb.txt:102 |
| 2 | `trace_width` | track_width.opt [netclass] | lossy: legacy board parser not in checkout (GK-K) | golden-rules/eez-psu-h24005.txt:37-38; golden-rules/esp32-evb.txt:102 |
| 3 | `via_dia` | via_diameter.opt [netclass] | lossy: legacy board parser not in checkout (GK-K) | golden-rules/eez-psu-h24005.txt:37-38; golden-rules/esp32-evb.txt:102 |
| 4 | `via_drill` | hole_size.opt [netclass, vias] | lossy: legacy board parser not in checkout (GK-K) | golden-rules/eez-psu-h24005.txt:37-38; golden-rules/esp32-evb.txt:102 |
| 5 | `uvia_dia` | via_diameter.opt [netclass, micro] | lossy: legacy board parser not in checkout (GK-K) | golden-rules/eez-psu-h24005.txt:37-38; golden-rules/esp32-evb.txt:102 |
| 6 | `uvia_drill` | hole_size.opt [netclass, micro] | lossy: legacy board parser not in checkout (GK-K) | golden-rules/eez-psu-h24005.txt:37-38; golden-rules/esp32-evb.txt:102 |
| 7 | `diff_pair_gap` | diff_pair_gap.opt [netclass] | lossy: legacy board parser not in checkout (GK-K) | golden-rules/eez-psu-h24005.txt:37-38; golden-rules/esp32-evb.txt:102 |
| 8 | `diff_pair_width` | track_width.opt [netclass, diffpair] | lossy: legacy board parser not in checkout (GK-K) | golden-rules/eez-psu-h24005.txt:37-38; golden-rules/esp32-evb.txt:102 |
| 9 | `add_net` | @netclass membership | lossy: legacy board parser not in checkout (GK-K) | golden-rules/eez-psu-h24005.txt:37-38; golden-rules/esp32-evb.txt:102 |

#### 4.1.17 KiCad 5 `.pro` `[pcbnew]` keys (seen in esp32-evb, glasgow, sh-esp32, speeduino-v04, eez-psu-h24005)

| # | .pro [pcbnew] key | Canonical target | Quality | Provenance |
|---|---|---|---|---|
| 1 | `CopperLayerCount` | stackup.layers | lossy: legacy .pro parser not in checkout; meaning unverified (GK-K) | golden/esp32-evb/src/HARDWARE/REV-F/ESP32-EVB_Rev_F.pro (golden-src) |
| 2 | `BoardThickness` | stackup.board_thickness | lossy: legacy .pro parser not in checkout; meaning unverified (GK-K) | golden/esp32-evb/src/HARDWARE/REV-F/ESP32-EVB_Rev_F.pro (golden-src) |
| 3 | `AllowMicroVias` | via_span.allowed_types (micro) | lossy: legacy .pro parser not in checkout; meaning unverified (GK-K) | golden/esp32-evb/src/HARDWARE/REV-F/ESP32-EVB_Rev_F.pro (golden-src) |
| 4 | `AllowBlindVias` | via_span.allowed_types (blind, buried) | lossy: legacy .pro parser not in checkout; meaning unverified (GK-K) | golden/esp32-evb/src/HARDWARE/REV-F/ESP32-EVB_Rev_F.pro (golden-src) |
| 5 | `RequireCourtyardDefinitions` | courtyard_required[missing] enabled | lossy: legacy .pro parser not in checkout; meaning unverified (GK-K) | golden/esp32-evb/src/HARDWARE/REV-F/ESP32-EVB_Rev_F.pro (golden-src) |
| 6 | `ProhibitOverlappingCourtyards` | courtyard_clearance enabled | lossy: legacy .pro parser not in checkout; meaning unverified (GK-K) | golden/esp32-evb/src/HARDWARE/REV-F/ESP32-EVB_Rev_F.pro (golden-src) |
| 7 | `MinTrackWidth` | track_width.min [board_floor] | lossy: legacy .pro parser not in checkout; meaning unverified (GK-K) | golden/esp32-evb/src/HARDWARE/REV-F/ESP32-EVB_Rev_F.pro (golden-src) |
| 8 | `MinViaDiameter` | via_diameter.min [board_floor] | lossy: legacy .pro parser not in checkout; meaning unverified (GK-K) | golden/esp32-evb/src/HARDWARE/REV-F/ESP32-EVB_Rev_F.pro (golden-src) |
| 9 | `MinViaDrill` | hole_size.min [board_floor, vias] | lossy: legacy .pro parser not in checkout; meaning unverified (GK-K) | golden/esp32-evb/src/HARDWARE/REV-F/ESP32-EVB_Rev_F.pro (golden-src) |
| 10 | `MinMicroViaDiameter` | via_diameter.min [board_floor, micro] | lossy: legacy .pro parser not in checkout; meaning unverified (GK-K) | golden/esp32-evb/src/HARDWARE/REV-F/ESP32-EVB_Rev_F.pro (golden-src) |
| 11 | `MinMicroViaDrill` | hole_size.min [board_floor, micro] | lossy: legacy .pro parser not in checkout; meaning unverified (GK-K) | golden/esp32-evb/src/HARDWARE/REV-F/ESP32-EVB_Rev_F.pro (golden-src) |
| 12 | `MinHoleToHole` | hole_to_hole.min [board_floor] | lossy: legacy .pro parser not in checkout; meaning unverified (GK-K) | golden/esp32-evb/src/HARDWARE/REV-F/ESP32-EVB_Rev_F.pro (golden-src) |
| 13 | `SolderMaskClearance` | solder_mask_expansion.opt [board default] | lossy: legacy .pro parser not in checkout; meaning unverified (GK-K) | golden/eez-psu-h24005/src/Consolidated/KiCad/EEZ PSU consolidated r5B13a.pro (golden-src) |
| 14 | `SolderMaskMinWidth` | solder_mask_bridge.min_web + solder_mask_sliver.min [board_floor] (as section 4.1.4 row 2) | lossy: legacy .pro parser not in checkout; meaning unverified (GK-K) | golden/eez-psu-h24005/src/Consolidated/KiCad/EEZ PSU consolidated r5B13a.pro (golden-src) |
| 15 | `SolderPasteClearance` | solder_paste_margin.abs [board default] | lossy: legacy .pro parser not in checkout; meaning unverified (GK-K) | golden/esp32-evb/src/HARDWARE/REV-F/ESP32-EVB_Rev_F.pro (golden-src) |
| 16 | `SolderPasteRatio` | solder_paste_margin.ratio [board default] | lossy: legacy .pro parser not in checkout; meaning unverified (GK-K) | golden/esp32-evb/src/HARDWARE/REV-F/ESP32-EVB_Rev_F.pro (golden-src) |
| 17 | `[pcbnew/Netclasses/..] sections` | netclass-origin rules (clearance, track_width, via_diameter, hole_size, diff_pair_gap) | lossy: legacy .pro parser not in checkout; meaning unverified (GK-K) | golden/esp32-evb/src/HARDWARE/REV-I/ESP32-EVB_Rev_I.pro (golden-src) |
| 18 | `TrackWidthN, ViaDiameterN, dPair*N` | router_presets | lossy: legacy .pro parser not in checkout; meaning unverified (GK-K) | golden/glasgow/src/hardware/boards/glasgow/revC1/glasgow.pro (golden-src) |
| 19 | `ERC_TestSimilarLabels` | erc_label[similar_labels] enabled | lossy: legacy .pro parser not in checkout; meaning unverified (GK-K) | golden/eez-psu-h24005/src/Consolidated/KiCad/EEZ PSU consolidated r5B13a.pro (golden-src) |

#### 4.1.18 Obsolete and legacy project keys seen in golden boards

| # | Key | Canonical target | Quality | Provenance |
|---|---|---|---|---|
| 1 | `rules.allow_microvias` | via_span.allowed_types (micro) | lossy: current code does not read it (Q-K9) | golden-rules/bms-c1.txt:3 |
| 2 | `rules.allow_blind_buried_vias` | via_span.allowed_types (blind, buried) | lossy: current code does not read it (Q-K9) | golden-rules/bms-c1.txt:3 |
| 3 | `rules.solder_mask_clearance` | solder_mask_expansion.opt [board default] | exact (migrated once to the board file) | golden-rules/esp32-evb.txt:3; pcbnew/board_design_settings.cpp:1150-1179 |
| 4 | `rules.solder_mask_min_width` | solder_mask_bridge.min_web + solder_mask_sliver.min [board_floor] (as section 4.1.4 row 2) | exact migration; lossy meaning (section 4.1.4 row 2) | golden-rules/esp32-evb.txt:3; pcbnew/board_design_settings.cpp:1150-1179 |
| 5 | `rules.solder_paste_clearance` | solder_paste_margin.abs [board default] | exact (migrated) | golden-rules/esp32-evb.txt:3; pcbnew/board_design_settings.cpp:1150-1179 |
| 6 | `rules.solder_paste_margin_ratio` | solder_paste_margin.ratio [board default] | exact (migrated) | golden-rules/esp32-evb.txt:3; pcbnew/board_design_settings.cpp:1150-1179 |
| 7 | `rule_severities.hole_near_hole` | @severity of hole_to_hole | exact (read before hole_to_hole) | golden-rules/glasgow.txt:4; pcbnew/board_design_settings.cpp:384-390 |
| 8 | `erc.rule_severities.global_label_dangling` | @severity of erc_label[single_global_label] (likely successor) | lossy: skipped by the current loader; successor unverified (Q-K9) | golden-rules/bms-c1.txt:15; eeschema/erc/erc_settings.cpp:155-169 |
| 9 | `tuning_profiles_impedance_geometric` | impedance_target profiles | lossy: ohm-to-geometry table; empty in the golden data | golden-rules/glasgow.txt:38 |

#### 4.1.19 Golden-board observations that affect the KiCad mapping

| # | Observation | Effect | Provenance |
|---|---|---|---|
| 1 | Projects still store `rules.allow_microvias` and `rules.allow_blind_buried_vias` | No reader in the checkout; the keys are ignored unless a project-file migration outside the checkout consumes them | golden-rules/bms-c1.txt:3; golden-rules/speeduino-v04.txt:3 |
| 2 | Projects store `rules.solder_mask_clearance`, `solder_mask_min_width`, `solder_paste_clearance`, `solder_paste_margin_ratio` | Read once by migration 1→2 into board-file members, then erased from the project | golden-rules/esp32-evb.txt:3; golden-rules/moteus.txt:3 |
| 3 | `pth_inside_courtyard` and `npth_inside_courtyard` set to `ignore` in 12 projects | Stored value wins. The current default is error, so a project without the key would now fail DRC on these | golden-rules/bms-c1.txt:4; golden-rules/hackrf-one.txt:4 |
| 4 | DRC severity key `hole_near_hole` | Read before `hole_to_hole`; a later `hole_to_hole` value overwrites it | golden-rules/glasgow.txt:4; pcbnew/board_design_settings.cpp:384-390 |
| 5 | ERC severity key `global_label_dangling` = warning in 15 projects | Not in the current ERC key list, so it is skipped on load; `ERCE_SINGLE_GLOBAL_LABEL` (key `single_global_label`) keeps its default, ignore. A migration may exist outside the checkout | golden-rules/bms-c1.txt:15; eeschema/erc/erc_settings.cpp:155-169, :105 |
| 6 | Every golden DRC severity key other than `hole_near_hole` matches a current DRC settings key | Keys round-trip | golden-rules/glasgow.txt:18; pcbnew/drc/drc_item.cpp:95-397 |
| 7 | Netclass JSON in golden projects uses exactly the keys in section 4.1.6; `tuning_profile` appears only in the 2026 Glasgow revisions | Keys round-trip | golden-rules/bms-c1.txt:11; golden-rules/glasgow.txt:25 |

#### 4.1.20 KiCad 6 and later board-file `(setup)` tokens (golden boards with file versions 20211014-20260206)

The 39 `(setup)` tokens that occur only in the newer golden board files. The 28 tokens these files share with KiCad 4
and 5 files (`pad_to_mask_clearance`, `aux_axis_origin`, `grid_origin`, `pcbplotparams` and its common children,
`viasonmask`, `subtractmaskfromsilk`, `drillshape`) resolve through section 4.1.15. Board-setup members are those of
section 4.1.4. Parser lines refer to the current s-expression loader, fetched at the same commit.

| # | (setup) token | Example | Unit | Canonical target | Quality | Note | Provenance |
|---|---|---|---|---|---|---|---|
| 1 | `allow_soldermask_bridges_in_footprints` | (key only) | bool | solder_mask_bridge.allow_in_footprint | exact | member of section 4.1.4 row 5 | golden-rules/bms-c1.txt:18; pcbnew/pcb_io/kicad_sexpr/pcb_io_kicad_sexpr_parser.cpp:3317-3320† |
| 2 | `tenting` | (key only) | flags (front, back) | via_protection.tent_top / tent_bottom [board default] | exact | section 4.1.4 row 6; written as bare `front back` symbols (older) or `(front yes) (back yes)` sub-lists (2026 files) | golden-rules/bms-c1.txt:18; golden-rules/glasgow.txt:229; pcbnew/pcb_io/kicad_sexpr/pcb_io_kicad_sexpr_parser.cpp:3322-3328† |
| 3 | `covering`, `plugging` | (key only) | flags (front, back) | via_protection.cover_top / cover_bottom, plug_top / plug_bottom [board default] | exact | section 4.1.4 row 7 | golden-rules/glasgow.txt:229; pcbnew/pcb_io/kicad_sexpr/pcb_io_kicad_sexpr_parser.cpp:3330-3344† |
| 4 | `capping`, `filling` | (key only) | bool | via_protection.cap, via_protection.fill [board default] | exact | section 4.1.4 rows 7-8 | golden-rules/glasgow.txt:229; pcbnew/pcb_io/kicad_sexpr/pcb_io_kicad_sexpr_parser.cpp:3346-3357† |
| 5 | `front`, `back` | (key only) | bool | the side flags of rows 2-3 | exact | children of `tenting`, `covering`, `plugging` | golden-rules/glasgow.txt:229; pcbnew/pcb_io/kicad_sexpr/pcb_io_kicad_sexpr_parser.cpp:3322-3344† |
| 6 | `pad_to_paste_clearance_ratio` | -0.12 | fraction | solder_paste_margin.ratio [board default] | exact | section 4.1.4 row 4 | golden-rules/hackrf-one.txt:14, :16; pcbnew/pcb_io/kicad_sexpr/pcb_io_kicad_sexpr_parser.cpp:3312-3315† |
| 7 | `stackup` | (block) | - | stackup | exact | section 3.11 | golden-rules/bms-c1.txt:18; include/board_design_settings.h:782; pcbnew/pcb_io/kicad_sexpr/pcb_io_kicad_sexpr_parser.cpp:3059-3064, :2226-2475† |
| 8 | `layer`, `type` (inside `stackup`) | (key only) | name, type string | stackup.layers: one entry per stackup layer with its name and type (copper, core, prepreg, silk, paste, mask) | exact | | golden-rules/bms-c1.txt:18; pcbnew/pcb_io/kicad_sexpr/pcb_io_kicad_sexpr_parser.cpp:2248, :2299, :2382-2384† |
| 9 | `thickness` (inside a stackup layer) | (key only) | mm | stackup per-layer thickness | exact | | golden-rules/bms-c1.txt:18; pcbnew/pcb_io/kicad_sexpr/pcb_io_kicad_sexpr_parser.cpp:2388-2399† |
| 10 | `material` (inside a stackup layer) | "JLCPCB Soldermask" | string | stackup per-layer material | exact | | golden-rules/glasgow.txt:230; pcbnew/pcb_io/kicad_sexpr/pcb_io_kicad_sexpr_parser.cpp:2406-2408† |
| 11 | `epsilon_r` (inside a stackup layer) | (key only) | dimensionless | stackup per-layer relative permittivity | exact | input to impedance_target geometry | golden-rules/bms-c1.txt:18; pcbnew/pcb_io/kicad_sexpr/pcb_io_kicad_sexpr_parser.cpp:2412-2414† |
| 12 | `loss_tangent` (inside a stackup layer) | (key only) | dimensionless | stackup per-layer loss tangent | exact | | golden-rules/bms-c1.txt:18; pcbnew/pcb_io/kicad_sexpr/pcb_io_kicad_sexpr_parser.cpp:2418-2420† |
| 13 | `color` (inside a stackup layer) | (key only) | string | stackup per-layer colour (solder mask or silkscreen colour, a fab order option) | exact | | golden-rules/glasgow.txt:229; pcbnew/pcb_io/kicad_sexpr/pcb_io_kicad_sexpr_parser.cpp:2451-2475† |
| 14 | `copper_finish` | (key only) | string | stackup.copper_finish (surface finish, a fab order option) | exact | | golden-rules/bms-c1.txt:18; pcbnew/pcb_io/kicad_sexpr/pcb_io_kicad_sexpr_parser.cpp:2252-2256† |
| 15 | `dielectric_constraints` | (key only) | bool | stackup.dielectric_constraints | exact | meaning inferred from the member name `m_HasDielectricConstrains` (the dielectric values are binding for the fab) | golden-rules/bms-c1.txt:18; pcbnew/pcb_io/kicad_sexpr/pcb_io_kicad_sexpr_parser.cpp:2264-2268† |
| 16 | `crossoutdnponfab`, `dashed_line_dash_ratio`, `dashed_line_gap_ratio`, `disableapertmacros`, `dxfimperialunits`, `dxfpolygonmode`, `dxfusepcbnewfont`, `hidednponfab`, `pdf_back_fp_property_popups`, `pdf_front_fp_property_popups`, `pdf_metadata`, `pdf_single_document`, `plot_black_and_white`, `plot_on_all_layers_selection`, `plotfptext`, `plotpadnumbers`, `sketchdnponfab`, `sketchpadsonfab`, `svgprecision`, `svguseinch` (20 children of `pcbplotparams`) | (key only) | - | none | none: plot setting, not a rule | | golden-rules/bms-c1.txt:18; golden-rules/hackrf-one.txt:14 (`svguseinch`); golden-rules/sh-esp32.txt:73 (`plotfptext`); common/pcb_plot_params.keywords:1-55†; pcbnew/pcb_plot_params.cpp:551-828† |

### 4.2 Altium -> canonical

#### 4.2.1 Rule kinds (all 70 `TRuleKind` ordinals 0-69)

Ordinals 0-51 and their strings come from Altium's own DelphiScript sample (W1:16-67); 32 kinds are confirmed by the uint16 prefix of each ODrive
`Rules6` record (golden/odrive-v3/src/v3/PCB.PcbDoc, 37 records, 19186 bytes consumed). Arity: `B`/`U` sourced,
`B*`/`U*` unverified (GK). Ordinals above 51 are unverified (L1: confirm via Altium scripting, Q-A12).

| # | Category | Rule type (UI) | RULEKIND (ordinal) | Arity | Canonical target | Quality | Note | Provenance |
|---|---|---|---|---|---|---|---|---|
| E1 | Electrical | Clearance | `Clearance` (0) | B | clearance.min (GAP); NETSCOPE -> pair; each OBJECTCLEARANCES entry -> extra clearance rule with type() pair (entries naming Hole -> hole_clearance); SameNetOnly -> physical_clearance with pair same_net; IGNOREPADTOPADCLEARANCEINFOOTPRINT -> ignore rule for pad pairs of one component; GENERICCLEARANCE (7.9921 mil, equal to GAP in ODrive) kept in params, meaning unverified (GK) | exact |  | params: golden-rules/odrive-v3.txt:35; W5:1; W6:1. Importer: pcbnew/pcb_io/altium/altium_parser_pcb.cpp:886-890 (ALTIUM_RULE_KIND::CLEARANCE, GAP only) |
| E2 | Electrical | Short-Circuit | `ShortCircuit` (15) | B* | short_circuit.allowed (+ tracks_crossing) | exact | Altium does not separate crossings from shorts | golden-rules/odrive-v3.txt:37; UI name AD4; arity unverified (GK) |
| E3 | Electrical | Un-Routed Net | `UnRoutedNet` (16) | U* | unrouted_net.check_bad_connections | exact |  | golden-rules/odrive-v3.txt:36; W6:27; UI name AD4 |
| E4 | Electrical | Un-Connected Pin | `UnConnectedPin` (45) | U* | unconnected_pin | lossy: keys and meaning unverified (GK) |  | UI name AD4; params and arity unverified (GK) |
| E5 | Electrical | Modified Polygon | `UnpouredPolygon` (62); enum name `eRule_ModifiedPolygon` (unverified (L1: confirm via Altium scripting, Q-A12)) | U* | unpoured_zone | exact | meaning unverified (GK) | golden-rules/odrive-v3.txt:8; UI name AD4; W13:95 |
| E6 | Electrical | Creepage Distance | `Creepage` (65) | B* | creepage | lossy: keys unknown |  | AD5; W13:180 (first and second scoped object, third-party, INFERRED); arity unverified (GK) |
| R1 | Routing | Width | `Width` (2) | U* | track_width.min/opt/max (MINLIMIT, PREFEREDWIDTH, MAXLIMIT; per-layer keys -> one rule per layer) + impedance_target (MINIMP, MAXIMP, FAVIMP, IMPEDANCEPROFILE*) | exact | CHECKCONNECTEDCOPPER kept in params (meaning unverified (GK)); mid-layer ordinals need the layer stack | params: golden-rules/odrive-v3.txt:31; golden/odrive-v3/src/v3/PCB.PcbDoc rec Width (MINIMP, MAXIMP, FAVIMP, CHECKCONNECTEDCOPPER, BOTTOMLAYER_PREFWIDTH); W5:2-6. Importer: pcbnew/pcb_io/altium/altium_parser_pcb.cpp:925-934 |
| R2 | Routing | Routing Topology | `RoutingTopology` (7) | U* | routing_topology.value | exact |  | golden-rules/odrive-v3.txt:24; arity unverified (GK) |
| R3 | Routing | Routing Priority | `RoutingPriority` (8) | U* | routing_priority.value | exact |  | golden-rules/odrive-v3.txt:30; arity unverified (GK) |
| R4 | Routing | Routing Layers | `RoutingLayers` (9) | U* | routing_layers.allowed_layers: one key per layer, `<LAYER NAME>_V5` = TRUE / FALSE (`TOP LAYER_V5`, `MID LAYER 1_V5` ... `MID LAYER 30_V5`, `BOTTOM LAYER_V5`); TRUE -> layer in allowed_layers | exact | layer names need the layer stack (Q-A11); the `_V5` suffix is kept verbatim in native | golden-rules/odrive-v3.txt:29; golden/odrive-v3/src/v3/PCB.PcbDoc rec RoutingLayers; arity unverified (GK) |
| R5 | Routing | Routing Corners | `RoutingCorners` (10) | U* | routing_corners.style, .setback | exact |  | golden-rules/odrive-v3.txt:28; arity unverified (GK) |
| R6 | Routing | Routing Via Style | `RoutingVias` (11) | U* | via_diameter.min/opt/max (MINWIDTH, WIDTH, MAXWIDTH) + hole_size.min/opt/max (MINHOLEWIDTH, HOLEWIDTH, MAXHOLEWIDTH) [scope type(via)] + via_span (VIASTYLE) | exact | via templates (USEVIATEMPLATES, VIATEMPLATE*) kept in params | params: golden-rules/odrive-v3.txt:27; W6:33. Importer: pcbnew/pcb_io/altium/altium_parser_pcb.cpp:915-924; pcbnew/pcb_io/altium/altium_pcb.cpp:3756-3760, :3791-3795, :186-189. Arity unverified (GK) |
| R7 | Routing | Fanout Control | `FanoutControl` (49) | U* | fanout_control | exact |  | golden-rules/odrive-v3.txt:16-19; W5:26-30 |
| R8 | Routing | Differential Pairs Routing | `DiffPairsRouting` (51) | U* | diff_pair_gap.min/opt/max (MINLIMIT, MOSTFREQGAP, MAXLIMIT) + diff_pair_uncoupled.max + track_width per layer [scope diffpair] + impedance_target | exact | gap reading of MIN/MAXLIMIT unverified (L1: confirm via Altium scripting, Q-A12) | params: golden-rules/odrive-v3.txt:15; golden/odrive-v3/src/v3/PCB.PcbDoc rec DiffPairsRouting (MAXUNCOUPLEDLENGTH, BOTTOMLAYER_*); W5:32-36. Importer: pcbnew/pcb_io/altium/altium_parser_pcb.cpp:896-899 (kind only, no parameters read; nothing in pcbnew/pcb_io/altium/*.cpp uses it) |
| R9 | Routing | Routing Neck-Down (name inferred) | `RoutingNeckDown` (67) | unknown | routing_neckdown | lossy: meaning unknown |  | Category and UI name inferred from the string only |
| S1 | SMT | SMD To Corner | `SMDToCorner` (23) | U* | smd_to_corner.min | lossy: unit of DISTANCE inferred (internal units) |  | W6:24; AD10; arity unverified (GK) |
| S2 | SMT | SMD To Plane | `SMDToPlane` (46) | U* | smd_to_plane.max | lossy: key unknown |  | AD10; arity unverified (GK) |
| S3 | SMT | SMD Neck-Down | `SMDNeckDown` (47) | U* | smd_neckdown.max | lossy: key unknown |  | AD10; arity unverified (GK) |
| S4 | SMT | SMD Entry | `SMDEntry` (60) | U* | smd_entry | exact |  | W6:4; AD10; arity unverified (GK) |
| M1 | Mask | Solder Mask Expansion | `SolderMaskExpansion` (13) | U* | solder_mask_expansion.opt (EXPANSION) + via_protection.tent_top (ISTENTINGTOP) | exact |  | params: golden-rules/odrive-v3.txt:25; W6:31; W5:16. Importer: pcbnew/pcb_io/altium/altium_parser_pcb.cpp:940-944 |
| M2 | Mask | Paste Mask Expansion | `PasteMaskExpansion` (14) | U* | solder_paste_margin.abs (EXPANSION) | exact |  | params: golden-rules/odrive-v3.txt:26. Importer: pcbnew/pcb_io/altium/altium_parser_pcb.cpp:935-939; pcbnew/pcb_io/altium/altium_pcb.cpp:3774, :3779-3780. UI name and arity unverified (GK) |
| P1 | Plane | Power Plane Connect Style | `PlaneConnect` (6) | U* | zone_connection (PLANECONNECTSTYLE) + thermal_relief_gap (RELIEFAIRGAP) + thermal_spoke_width (RELIEFCONDUCTORWIDTH) + thermal_spoke_count.count (RELIEFENTRIES) [scope type(zone.plane)] | lossy: RELIEFEXPANSION meaning unverified (GK), kept in params |  | params: golden-rules/odrive-v3.txt:34; W5:9 |
| P2 | Plane | Power Plane Clearance | `PlaneClearance` (12) | U* | plane_clearance.min | exact |  | params: golden-rules/odrive-v3.txt:33. Importer: pcbnew/pcb_io/altium/altium_parser_pcb.cpp:945-949; pcbnew/pcb_io/altium/altium_pcb.cpp:3612, :3618-3636. Arity unverified (GK) |
| P3 | Plane | Polygon Connect Style | `PolygonConnect` (20) | B* | zone_connection (CONNECTSTYLE) + thermal_relief_gap (AIRGAPWIDTH) + thermal_spoke_width (RELIEFCONDUCTORWIDTH) + thermal_spoke_count (RELIEFENTRIES, POLYGONRELIEFANGLE) | exact | THPAD./SMDPAD./VIA. prefixes become type() scopes | params: golden-rules/odrive-v3.txt:32; W5:20. Importer: pcbnew/pcb_io/altium/altium_parser_pcb.cpp:950-964; pcbnew/pcb_io/altium/altium_pcb.cpp:3638-3665 (TODO at :3658) |
| T1 | Testpoint | Fabrication Testpoint Style | `FabricationTestpoint` (43); `Testpoint` in W1:59 | U* | testpoint_style[fabrication]: MINSIZE / PREFEREDSIZE / MAXSIZE -> size min / opt / max; MINHOLESIZE / PREFEREDHOLESIZE / MAXHOLESIZE -> hole min / opt / max; TESTPOINTGRID -> grid; GRIDTOLERANCE -> grid_tolerance; USEGRID -> options.use_grid; ALLOWSIDETOP / ALLOWSIDEBOTTOM -> sides; TESTPOINTUNDERCOMPONENT -> under_component; SIDE kept in params | lossy: SIDE=3 unresolved; key meanings read from the names (GK) |  | golden-rules/odrive-v3.txt:4; string change: W1:59 vs unverified (L1: confirm via Altium scripting, Q-A12) and the ODrive record |
| T2 | Testpoint | Fabrication Testpoint Usage | `FabricationTestPointUsage` (44); `TestPointUsage` in W1:60 | U* | testpoint_usage[fabrication] | lossy: VALID codes from enum order only |  | golden-rules/odrive-v3.txt:7 |
| T3 | Testpoint | Assembly Testpoint Style | `AssemblyTestpoint` (57) | U* | testpoint_style[assembly]: the keys of T1 except SIDE (MINSIZE, PREFEREDSIZE, MAXSIZE, MINHOLESIZE, PREFEREDHOLESIZE, MAXHOLESIZE, TESTPOINTGRID, GRIDTOLERANCE, USEGRID, ALLOWSIDETOP, ALLOWSIDEBOTTOM, TESTPOINTUNDERCOMPONENT), mapped as in T1 | exact | key meanings read from the names (GK) | golden-rules/odrive-v3.txt:5; W6:6 |
| T4 | Testpoint | Assembly Testpoint Usage | `AssemblyTestPointUsage` (58) | U* | testpoint_usage[assembly] | exact |  | golden-rules/odrive-v3.txt:6 |
| F1 | Manufacturing | Minimum Annular Ring | `MinimumAnnularRing` (19) | U* | annular_width.min (MINIMUMRING) | exact |  | golden-rules/odrive-v3.txt:3; W6:26 |
| F2 | Manufacturing | Acute Angle | `AcuteAngle` (21) | U* | track_angle.min | lossy: keys unknown |  | UI name AD2; params unverified (GK) |
| F3 | Manufacturing | Hole Size | `HoleSize` (42) | U* | hole_size.min/max (absolute) or .min_pct/.max_pct (relative) | exact | percent reading unverified (GK) | golden-rules/odrive-v3.txt:20; W13:171. Importer: pcbnew/pcb_io/altium/altium_parser_pcb.cpp:904-909; pcbnew/pcb_io/altium/altium_pcb.cpp:3762-3765 |
| F4 | Manufacturing | Layer Pairs | `LayerPairs` (48) | U* | via_span.enforce_pairs (ENFORCE) | exact | meaning unverified (GK) | golden-rules/odrive-v3.txt:21 |
| F5 | Manufacturing | Hole To Hole Clearance | `HoleToHoleClearance` (52) | B* | hole_to_hole.min (GAP) + microvia_stack.allow_stacked (ALLOWSTACKEDMICROVIAS) | exact |  | golden-rules/odrive-v3.txt:14. Importer: pcbnew/pcb_io/altium/altium_parser_pcb.cpp:910-914; pcbnew/pcb_io/altium/altium_pcb.cpp:3744, :3767-3768 |
| F6 | Manufacturing | Minimum Solder Mask Sliver | `MinimumSolderMaskSliver` (53) | B* | solder_mask_sliver.min (MINSOLDERMASKWIDTH) | exact |  | golden-rules/odrive-v3.txt:13 |
| F7 | Manufacturing | Silk To Solder Mask Clearance | `SilkToSolderMaskClearance` (54) | B* | silk_to_mask_clearance.min, .measure_to (CLEARANCETOEXPOSEDCOPPER) | exact |  | golden-rules/odrive-v3.txt:12; W13:178 |
| F8 | Manufacturing | Silk To Silk Clearance | `SilkToSilkClearance` (55) | B* | silk_clearance.min (SILKTOSILKCLEARANCE) | exact |  | golden-rules/odrive-v3.txt:11 |
| F9 | Manufacturing | Net Antennae | `NetAntennae` (56) | U* | dangling_track.max_stub (NETANTENNAETOLERANCE) | exact |  | golden-rules/odrive-v3.txt:9 |
| F10 | Manufacturing | Silk To Board Region Clearance | `SilkToBoardRegionClearance` (59) | B* | silk_to_edge_clearance | lossy: value key unknown |  | golden-rules/odrive-v3.txt:10; W6:41; category unverified (GK) |
| F11 | Manufacturing | Board Outline Clearance | `BoardOutlineClearance` (63) | B* | edge_clearance.min (GAP; OBJECTCLEARANCES -> extra rules) | exact |  | W6:2; W11:1; W7 rec BoardOutlineClearance. Importer: pcbnew/pcb_io/altium/altium_parser_pcb.cpp:891-895 |
| H1 | High Speed | Parallel Segment | `ParallelSegment` (1) | B* | parallel_segment | lossy: keys unknown |  | AD8; arity unverified (GK) |
| H2 | High Speed | Length | `Length` (3) | U* | length.min/max | lossy: keys unknown |  | UI name AD4; params unverified (GK) |
| H3 | High Speed | Matched Lengths | `MatchedLengths` (4) | U* | skew.max (TOLERANCE) + skew.compare (CHECK* flags) | exact |  | W5:7-8; AD8; arity unverified (GK) |
| H4 | High Speed | Daisy Chain Stub Length | `StubLength` (5); `eRule_DaisyChainStubLength` | U* | stub_length.max | lossy: key unknown; Altium daisy-chain stub differs from KiCad chain stub |  | AD8; arity unverified (GK) |
| H5 | High Speed | Vias Under SMD | `ViasUnderSMD` (17) | U* | vias_under_smd.allowed | exact |  | W7 rec ViasUnderSMD; UI name AD4 |
| H6 | High Speed | Maximum Via Count | `MaximumViaCount` (18) | U* | via_count.max | lossy: key unknown |  | UI name AD4 (unverified key) |
| H7 | High Speed | Max Via Stub Length (Back Drilling) | `BackDrilling` (64) | U* | backdrill_stub | lossy: keys unknown |  | AD7; arity unverified (GK) |
| H8 | High Speed | Return Path | `ReturnPath` (66) | U | return_path | lossy: keys unknown |  | AD6 (arity "unary-classified" is from AD6) |
| L1 | Placement | Room Definition | `RoomDefinition` (22); `eRule_ConfinementConstraint` | U* | room_confinement (region, mode, lock) | exact |  | W8 rec RoomDefinition; W10:1; AD11 |
| L2 | Placement | Component Clearance | `ComponentClearance` (24) | B* | courtyard_clearance.min (GAP), .vertical_min (VERTICALGAP), .check_mode (COLLISIONCHECKMODE) | lossy: 3D bodies vs 2D courtyards; COLLISIONCHECKMODE=3 unresolved |  | golden-rules/odrive-v3.txt:23, :38, :39; W13:179 |
| L3 | Placement | Component Orientations | `ComponentOrientations` (25); `eRule_ComponentRotations` | U* | component_orientation | lossy: keys unknown |  | UI name AD4; params unverified (GK) |
| L4 | Placement | Permitted Layers | `PermittedLayers` (26) | U* | permitted_layers | lossy: keys unknown |  | UI name AD4; params unverified (GK) |
| L5 | Placement | Nets to Ignore | `NetsToIgnore` (27) | U* | nets_to_ignore | lossy: meaning unverified (GK) |  | UI name AD4; params unverified (GK) |
| L6 | Placement | Height | `Height` (50) | U* | component_height.min/opt/max | exact |  | golden-rules/odrive-v3.txt:22; AD11. Importer: pcbnew/pcb_io/altium/altium_parser_pcb.cpp:900-903 |
| I1 | Signal Integrity | Signal Stimulus | `SignalStimulus` (28) | U* | si_stimulus | lossy: keys unknown |  | AD9 |
| I2 | Signal Integrity | Overshoot - Falling Edge | `OvershootFalling` (29) | U* | si_overshoot[falling] | lossy: keys unknown |  | AD9 |
| I3 | Signal Integrity | Overshoot - Rising Edge | `OvershootRising` (30) | U* | si_overshoot[rising] | lossy: keys unknown |  | AD9 |
| I4 | Signal Integrity | Undershoot - Falling Edge | `UndershootFalling` (31) | U* | si_undershoot[falling] | lossy: keys unknown |  | AD9 |
| I5 | Signal Integrity | Undershoot - Rising Edge | `UndershootRising` (32) | U* | si_undershoot[rising] | lossy: keys unknown |  | AD9 |
| I6 | Signal Integrity | Impedance | `MaxMinImpedance` (33) | U* | si_impedance | lossy: keys unknown; unit stated as V in AD9 (likely ohm) |  | AD9 |
| I7 | Signal Integrity | Signal Top Value | `SignalTopValue` (34) | U* | si_signal_top_value | lossy: keys unknown |  | AD9 |
| I8 | Signal Integrity | Signal Base Value | `SignalBaseValue` (35) | U* | si_signal_base_value | lossy: keys unknown |  | AD9 |
| I9 | Signal Integrity | Flight Time - Rising Edge | `FlightTimeRising` (36) | U* | si_flight_time[rising] | lossy: keys unknown |  | AD9 |
| I10 | Signal Integrity | Flight Time - Falling Edge | `FlightTimeFalling` (37) | U* | si_flight_time[falling] | lossy: keys unknown |  | AD9 |
| I11 | Signal Integrity | Slope - Rising Edge | `SlopeRising` (39) | U* | si_slope[rising] | lossy: keys unknown |  | AD9 |
| I12 | Signal Integrity | Slope - Falling Edge | `SlopeFalling` (40) | U* | si_slope[falling] | lossy: keys unknown |  | AD9 |
| I13 | Signal Integrity | Supply Nets | `SupplyNets` (41) | U* | supply_net_voltage.value (VOLTAGE) | exact |  | W9 rec SupplyNets; AD9; W13:143 (generation from power ports, third-party) |
| X1 | outside categories | Layer Stack (legacy) | `LayerStack` (38) | unknown | stackup | lossy: keys unknown |  | W1:54; absence from category lists per AD4 summary |
| X2 | outside categories | none | no string (61); `eRule_None` | n/a | none | none: ordinal 61 (eRule_None) is a placeholder, not a rule |  | odrive-DRC |
| X3 | outside categories | Wire Bonding (name inferred) | `WireBonding` (68) | unknown | wire_bonding | lossy: meaning and keys unknown |  | AD12 |
| X4 | outside categories | Z-Axis Clearance | `ZAxisClearance` (69) | B* | z_axis_clearance.min | lossy: keys unknown |  | AD12 |

#### 4.2.2 Rule record fields

| # | Field | Canonical | Quality | Note | Provenance |
|---|---|---|---|---|---|
| C1 | `RULEKIND` | @kind (via the kind table) + @native.kind | model |  | pcbnew/pcb_io/altium/altium_parser_pcb.cpp:885-967; pcbnew/pcb_io/altium/altium_parser_pcb.h:111-127 |
| C2 | `NAME` | @native.name (+ @id) | model |  | pcbnew/pcb_io/altium/altium_parser_pcb.cpp:878; golden-rules/odrive-v3.txt:17, :38 |
| C3 | `PRIORITY` | @priority (per kind; 1 = highest) | model |  | pcbnew/pcb_io/altium/altium_parser_pcb.cpp:879; pcbnew/pcb_io/altium/altium_pcb.cpp:3729-3738 |
| C4 | `ENABLED` | @enabled | model |  | pcbnew/pcb_io/altium/altium_parser_pcb.cpp:880; pcbnew/pcb_io/altium/altium_pcb.cpp:1546, :1562, :161 |
| C5 | `SCOPE1EXPRESSION`, `SCOPE2EXPRESSION` | @scope.expr (SCOPE1), @scope.expr_b (SCOPE2) | model | section 2.5 | pcbnew/pcb_io/altium/altium_parser_pcb.cpp:882-883; pcbnew/pcb_io/altium/altium_parser_pcb.h:658-659 |
| C6 | `NETSCOPE` | @scope.pair: DifferentNets -> different_net, SameNetOnly -> same_net, AnyNet -> none | model |  | golden-rules/odrive-v3.txt:35; W6:1 |
| C7 | `LAYERKIND` | @scope.pair: SameLayer -> same_layer, AdjacentLayer -> adjacent_layer | model |  | all files (enum SameLayer, AdjacentLayer) |
| C8 | `COMMENT` | @description | model |  | W6:31 |
| C9 | `UNIQUEID` | @native (unique id) | model |  | golden/odrive-v3/src/v3/PCB.PcbDoc rec Clearance |
| C10 | `DEFINEDBYLOGICALDOCUMENT` | @tags+=defined_by_schematic | model |  | W9 rec SupplyNets; FALSE in all 37 ODrive records (`golden/odrive-v3/src/v3/PCB.PcbDoc rec`, every NAME) |
| C11 | `SELECTION`, `LAYER`, `LOCKED`, `POLYGONOUTLINE`, `USERROUTED`, `KEEPOUT`, `UNIONINDEX` | none | none: generic primitive fields with no rule meaning |  | golden/odrive-v3/src/v3/PCB.PcbDoc rec MinimumAnnularRing; W9 rec SupplyNets |
| C12 | Legacy `SCOPE<i>COUNT`, `SCOPE<i>_<j>_KIND`, `SCOPE<i>_<j>_VALUE`, `COMMENTLENGTH` | @scope.expr (legacy scope list, section 2.5) | model | the exporter writes both forms | W7 rec BoardOutlineClearance; W8 rec RoomDefinition |

#### 4.2.3 Scope queries

See section 2.5 (query terms, operators, legacy scope kinds, object ids). Priority: `PRIORITY` is an integer per
kind, 1 highest; the applied rule is the highest-priority enabled rule whose scopes match (AD13;
pcbnew/pcb_io/altium/altium_parser_pcb.h:1027-1031). Section 5 has the precedence model.

#### 4.2.4 PCB DRC options stream (`Design Rule Checker Options6/Data`)

D1 and D2 list `TRuleKind` ordinals that batch and online DRC check (ODrive: batch 0, 1, 2, 3, 4, 5, 6, 15, 16, 18,
19, 21, 22, 23, 42, 45, 46, 47, 50, 52, 53, 54, 55, 56, 59, 60, 62; online 0, 1, 2, 3, 4, 5, 9, 11, 15, 17, 18, 22,
23, 24, 45, 46, 47, 50, 51, 55, 60). They become `enabled` per kind for the two DRC modes.

| # | Key | Value (ODrive) | Canonical | Quality | Provenance |
|---|---|---|---|---|---|
| D1 | `RULESETTOCHECK` (batch DRC) | 0,1,2,3,4,5,6,15,16,18,19,21,22,23,42,45,46,47,50,52,53,54,5 | @enabled per kind (batch DRC) | model | odrive-DRC D1 |
| D2 | `ONLINERULESETTOCHECK` | 0,1,2,3,4,5,9,11,15,17,18,22,23,24,45,46,47,50,51,55,60 | @enabled per kind (online DRC) | model | odrive-DRC D2 |
| D3 | `MAXVIOLATIONCOUNT` | 50000 | drc_engine_settings.max_violations | exact | odrive-DRC D3 |
| D4 | `REPORTDRILLEDSMTPADS`, `REPORTINVALIDMULTILAYERPADS`, `INTERNALPLANEWARNINGS`, `VERIFYSHORTINGCOPPER`, `REPORTBROKENPLANES` | TRUE | footprint_integrity[drilled_smd_pad, invalid_multilayer_pad] + plane_integrity[internal_plane_warning] + short_circuit (VERIFYSHORTINGCOPPER) + plane_integrity[broken_plane] | lossy: options, meanings unverified (GK) | odrive-DRC D4 |
| D5 | `REPORTDEADCOPPER`, `DEADCOPPERMINAREA` | TRUE, 10000000000.000000 | isolated_copper.min_area | lossy: area unit unknown | odrive-DRC D5 |
| D6 | `REPORTSTARVEDTHERMALS`, `MINSTARVEDCOPPERPERCENT` | TRUE, 50 | starved_thermal.min_copper (%) | lossy: DRC option, not a rule; meaning of MINSTARVEDCOPPERPERCENT unverified (GK), only the value 50 is on disk | odrive-DRC D6 |
| D7 | `REPORTSTRADLINGHOLES`, `REPORTHOLESINVOIDS` | FALSE | plane_integrity[straddling_hole, hole_in_void] | lossy: meanings unverified (GK) | odrive-DRC D7 |
| D8 | `DOMAKEDRCFILE`, `DOMAKEDRCERRORLIST`, `DOSUBNETDETAILS`, `CHECKEXTERNALNETLIST`, `EXTERNALNETLISTFILENAME`, `REPORTFILENAME` | FALSE, TRUE, TRUE, FALSE, empty, empty | drc_engine_settings (report output) | model | odrive-DRC D8 |

#### 4.2.5 Schematic side (`.PrjPcb`)

| # | Item | Canonical | Quality | Provenance |
|---|---|---|---|---|
| N1 | Where stored | ERC rule set (stored in the project file) | model | golden/odrive-v3/src/v3/two_ax_PCB.PrjPcb:4507, :4620 |
| N2 | Error reporting | @severity of erc_undriven_input, erc_unconnected_pin, erc_no_connect, erc_dangling, erc_label, erc_bus, erc_annotation, erc_multi_unit, erc_symbol_library, erc_power, erc_netclass, erc_drawing, erc_simulation_model | lossy: Type&lt;k> to check map unknown (Q-A6) | golden/odrive-v3/src/v3/two_ax_PCB.PrjPcb:4508-4618; AD14 |
| N3 | Connection matrix | erc_pin_type_conflict.matrix | lossy: row and column labels unverified (GK) | golden/odrive-v3/src/v3/two_ax_PCB.PrjPcb:4621-4637; AD14 |
| N4 | How it is applied | ERC semantics (section 5) | model | AD14 |
| N5 | Related sections, not ERC | board_schematic_parity (Difference Levels); @netclass and component-class generation | lossy: meanings unverified (GK) | golden/odrive-v3/src/v3/two_ax_PCB.PrjPcb:4378, :4460, :4650-4657 |
| N6 | Schematic-defined PCB rules | @tags+=defined_by_schematic | model | W9 rec SupplyNets |

#### 4.2.6 What KiCad's Altium importer keeps, and what canonical import keeps instead

KiCad's importer keeps 9 of the 70 kinds in some form (HoleSize is read but its values are discarded,
pcbnew/pcb_io/altium/altium_pcb.cpp:3743, :3762-3765), only for `All`/`All`, `InNetClass('x')`/`All` and polygon
scopes, never translates queries and never writes `.kicad_dru` rules (pcbnew/pcb_io/altium/altium_pcb.cpp:3729-3819;
pcbnew/pcb_io/altium/altium_parser_pcb.cpp:885-967). It is not a semantic reference.

KiCad's importer turns the `RULEKIND` string into one of 13 `ALTIUM_RULE_KIND` values (enum at
pcbnew/pcb_io/altium/altium_parser_pcb.h:111-127; string tests at pcbnew/pcb_io/altium/altium_parser_pcb.cpp:886-966). The table ties each value to its
row in section 4.2.1 and says whether the importer uses it:

| # | `ALTIUM_RULE_KIND` | RULEKIND string | Section 4.2.1 row | Parameters KiCad reads | KiCad importer use | Canonical import | Provenance |
|---|---|---|---|---|---|---|---|
| A0 | `UNKNOWN` (0) | every other string: 57 kinds (69 of the 70 ordinals have a string, ordinal 61 has none), e.g. `ShortCircuit`, `RoutingLayers`, `FanoutControl` | all rows not listed below | none | dropped | kept: mapped through section 4.2.1 like any other kind | pcbnew/pcb_io/altium/altium_parser_pcb.cpp:964-966 |
| A1 | `CLEARANCE` (1) | `Clearance` | E1 | GAP | K1, K7, K8 | E1 | pcbnew/pcb_io/altium/altium_parser_pcb.cpp:886-890 |
| A2 | `DIFF_PAIR_ROUTINGS` (2) | `DiffPairsRouting` | R8 | none | parsed, then dropped: no use in pcbnew/pcb_io/altium/*.cpp | R8 | pcbnew/pcb_io/altium/altium_parser_pcb.cpp:896-899 |
| A3 | `HEIGHT` (3) | `Height` | L6 | none | parsed, then dropped: no use in pcbnew/pcb_io/altium/*.cpp | L6 | pcbnew/pcb_io/altium/altium_parser_pcb.cpp:900-903 |
| A4 | `HOLE_SIZE` (4) | `HoleSize` | F3 | MINLIMIT, MAXLIMIT | read as the default rule, values discarded (TODO) | F3 | pcbnew/pcb_io/altium/altium_parser_pcb.cpp:904-909; pcbnew/pcb_io/altium/altium_pcb.cpp:3743, :3762-3765 |
| A5 | `HOLE_TO_HOLE_CLEARANCE` (5) | `HoleToHoleClearance` | F5 | GAP | K4 | F5 | pcbnew/pcb_io/altium/altium_parser_pcb.cpp:910-914 |
| A6 | `WIDTH` (6) | `Width` | R1 | MINLIMIT, MAXLIMIT, PREFEREDWIDTH | K2, K7 | R1 | pcbnew/pcb_io/altium/altium_parser_pcb.cpp:925-934 |
| A7 | `PASTE_MASK_EXPANSION` (7) | `PasteMaskExpansion` | M2 | EXPANSION | K6 | M2 | pcbnew/pcb_io/altium/altium_parser_pcb.cpp:935-939 |
| A8 | `SOLDER_MASK_EXPANSION` (8) | `SolderMaskExpansion` | M1 | EXPANSION | K6 | M1 | pcbnew/pcb_io/altium/altium_parser_pcb.cpp:940-944 |
| A9 | `PLANE_CLEARANCE` (9) | `PlaneClearance` | P2 | CLEARANCE | K8 | P2 | pcbnew/pcb_io/altium/altium_parser_pcb.cpp:945-949 |
| A10 | `POLYGON_CONNECT` (10) | `PolygonConnect` | P3 | AIRGAPWIDTH, RELIEFCONDUCTORWIDTH, RELIEFENTRIES, CONNECTSTYLE | K9 | P3 | pcbnew/pcb_io/altium/altium_parser_pcb.cpp:950-963 |
| A11 | `ROUTING_VIAS` (11) | `RoutingVias` | R6 | WIDTH, MINWIDTH, MAXWIDTH, HOLEWIDTH, MINHOLEWIDTH, MAXHOLEWIDTH | K3, K7 | R6 | pcbnew/pcb_io/altium/altium_parser_pcb.cpp:915-924 |
| A12 | `BOARD_OUTLINE_CLEARANCE` (12) | `BoardOutlineClearance` | F11 | GAP | K5 | F11 | pcbnew/pcb_io/altium/altium_parser_pcb.cpp:891-895 |

| # | Altium source | KiCad importer target | Canonical import instead | Provenance |
|---|---|---|---|---|
| K1 | Clearance default `GAP` | Board minimum clearance; Default netclass clearance | canonical keeps the All/All Clearance as an explicit rule; KiCad folds it into board floor and Default netclass | pcbnew/pcb_io/altium/altium_pcb.cpp:3747-3748, :3785-3786 |
| K2 | Width default `MINLIMIT` / `PREFEREDWIDTH` | Board minimum track width / Default netclass track width | canonical keeps MINLIMIT, PREFEREDWIDTH and MAXLIMIT; KiCad drops MAXLIMIT | pcbnew/pcb_io/altium/altium_pcb.cpp:3750-3754, :3788-3789 |
| K3 | Routing Vias default `MINWIDTH`, `MINHOLEWIDTH` / `WIDTH`, `HOLEWIDTH` | Board minimum via diameter and through-hole drill / Default netclass via diameter and drill | canonical keeps min, opt and max of diameter and drill; KiCad drops the maxima | pcbnew/pcb_io/altium/altium_pcb.cpp:3756-3760, :3791-3795 |
| K4 | Hole To Hole default `GAP` | Board minimum hole-to-hole | canonical keeps scoped rules and ALLOWSTACKEDMICROVIAS; KiCad keeps the All/All GAP only | pcbnew/pcb_io/altium/altium_pcb.cpp:3767-3768 |
| K5 | Board Outline default `GAP` | Board copper-to-edge clearance | canonical keeps OBJECTCLEARANCES entries; KiCad keeps the All/All GAP only | pcbnew/pcb_io/altium/altium_pcb.cpp:3770-3771 |
| K6 | Solder / Paste Mask default `EXPANSION` | Board solder mask expansion / solder paste margin | canonical keeps scoped mask and paste rules (for example the IsVia tenting rule); KiCad keeps the All/All value only | pcbnew/pcb_io/altium/altium_pcb.cpp:3773-3780 |
| K7 | Netclass-scoped Clearance, Width, Routing Vias | That netclass's clearance, track width, via diameter and drill | canonical keeps any query scope; KiCad handles InNetClass('x')/All only | pcbnew/pcb_io/altium/altium_pcb.cpp:142-195 |
| K8 | Polygon-chosen Plane Clearance and Clearance | Zone local clearance: the larger of the two, Plane Clearance only if the zone covers a power or mixed layer, Clearance only if it covers a signal or mixed layer | canonical keeps PlaneClearance as plane_clearance; KiCad turns it into a zone local clearance | pcbnew/pcb_io/altium/altium_pcb.cpp:3612-3636 |
| K9 | Polygon-chosen Polygon Connect | Zone pad connection, spoke width, thermal gap, minimum thickness clamp | canonical keeps RELIEFENTRIES and the angle; KiCad parses but drops RELIEFENTRIES | pcbnew/pcb_io/altium/altium_pcb.cpp:3638-3665 |

### 4.3 Eagle -> canonical

#### 4.3.1 Where Eagle keeps rule data

| # | Container | Canonical | Quality | Note | Provenance |
|---|---|---|---|---|---|
| 1 | `<board><designrules name="...">` in a `.brd` | rule set for the board (all &lt;param> rows) | model |  | pcbnew/pcb_io/eagle/pcb_io_eagle.cpp:544, :590, :620-632, :162-229 |
| 2 | Standalone `.dru` file | same rule set from a standalone file | model | KiCad ignores .dru files; canonical import reads them | golden-rules/moteus.txt:114-115; pcbnew/pcb_io/eagle/pcb_io_eagle.cpp:521-617, :620-632 |
| 3 | `<board><classes><class>` | @netclass definitions + class matrix rules | model |  | common/io/eagle/eagle_parser.h:1487-1507; common/io/eagle/eagle_parser.cpp:1870-1901 |
| 4 | `<signal class="N">` | @netclass membership | model |  | pcbnew/pcb_io/eagle/pcb_io_eagle.cpp:2779-2789 |
| 5 | Per-object attributes | item_override rules and object data (table 4.3.4) | model |  | common/io/eagle/eagle_parser.h:961-1004, :550-575, :854-855, :876, :911-922 |
| 6 | Restrict layers tRestrict, bRestrict, vRestrict (Eagle layers 41, 42, 43) | disallow [scope in_area(restrict shape)] | exact |  | common/io/eagle/eagle_parser.h:1240-1242; pcbnew/pcb_io/eagle/pcb_io_eagle.cpp:121-159 |
| 7 | Eagle 4/5 binary board: a fixed 244-byte DRC block after the notes section | same parameters (binary DRC block) | lossy: only 5 values are recoverable |  | common/io/eagle/eagle_bin_parser.cpp:1286-1346, :2107-2121, :2543-2551 |
| 8 | Schematic `<classes>` | @netclass definitions (schematic side) | lossy: no consumer in the checkout (Q-E7) |  | common/io/eagle/eagle_parser.h:2046; common/io/eagle/eagle_parser.cpp:2723-2735 |

#### 4.3.2 `<param>` names (all 74 rows: 72 parameters plus the rule-set name and description)

| # | Parameter | Unit / form | Canonical target | Quality | Note | Provenance |
|---|---|---|---|---|---|---|
| 1 | `name` attribute of `<designrules>` | string | @ruleset.name (@native.name) | model | a trailing ' *' is kept verbatim; its meaning is a data hypothesis | golden-rules/eez-psu-h24005.txt:47; golden-rules/moteus.txt:90, :102, :106 |
| 2 | description (`<description language>` in `.brd`; `description[en]` in `.dru`) | text | @ruleset.description | model |  | golden-rules/moteus.txt:115; [X] raw :12194-12200 |
| 3 | `layerSetup` | string over layer numbers 1-16 | stackup.layers + via_span.allowed_pairs | lossy: syntax unverified (GK-E) |  | golden-rules/eez-psu-h24005.txt:48; golden-rules/moteus.txt:91, :115 |
| 4 | `mtCopper` | 16 values, mm | stackup copper thickness per layer | lossy: indexing hypothesis (GK-E) |  | golden-rules/eez-psu-h24005.txt:48; golden-rules/moteus.txt:91 |
| 5 | `mtIsolate` | 15 values, mm | stackup dielectric thickness per gap | lossy: indexing hypothesis (GK-E) |  | golden-rules/eez-psu-h24005.txt:48; golden-rules/moteus.txt:91 |
| 6 | `mdWireWire` | length | clearance.min [board_floor, type(track) vs type(track)] | exact | meaning from OMX | golden-rules/eez-psu-h24005.txt:48; golden-rules/moteus.txt:91, :103, :107, :115 |
| 7 | `mdWirePad` | length | clearance.min [board_floor, type(track) vs type(pad.pth)] | lossy: pad attribute split assumed; meaning unverified (GK-E) |  | golden-rules/eez-psu-h24005.txt:48; golden-rules/moteus.txt:91, :115 |
| 8 | `mdWireVia` | length | clearance.min [board_floor, type(track) vs type(via)] | exact | meaning from OMX | as row 7; [X] raw :12206. Meaning: OMX, as row 7. Import: pcbnew/pcb_io/eagle/pcb_io_eagle.cpp:175-224 |
| 9 | `mdPadPad` | length | clearance.min [board_floor, type(pad.pth) vs type(pad.pth)] | lossy: meaning unverified (GK-E) |  | as row 7; [X] raw :12207. Meaning: OMX, as row 7. Import: pcbnew/pcb_io/eagle/pcb_io_eagle.cpp:175-224 |
| 10 | `mdPadVia` | length | clearance.min [board_floor, type(pad.pth) vs type(via)] | lossy: meaning unverified (GK-E) |  | as row 7; [X] raw :12208. Meaning: OMX, as row 7. Import: pcbnew/pcb_io/eagle/pcb_io_eagle.cpp:175-224 |
| 11 | `mdViaVia` | length | clearance.min [board_floor, type(via) vs type(via)] | exact | meaning from OMX | as row 7; [X] raw :12209. Meaning: OMX, as row 7. Import: pcbnew/pcb_io/eagle/pcb_io_eagle.cpp:175-224 |
| 12 | `mdSmdPad` | length | clearance.min [board_floor, type(pad.smd) vs type(pad.pth)] | lossy: 0 mil may mean unchecked; meaning unverified (GK-E) |  | golden-rules/eez-psu-h24005.txt:48; golden-rules/moteus.txt:91 |
| 13 | `mdSmdVia` | length | clearance.min [board_floor, type(pad.smd) vs type(via)] | lossy: meaning unverified (GK-E) |  | as row 12; [X] raw :12211. Meaning: OMX, as row 7. Import: pcbnew/pcb_io/eagle/pcb_io_eagle.cpp:175-224 |
| 14 | `mdSmdSmd` | length | clearance.min [board_floor, type(pad.smd) vs type(pad.smd)] | lossy: meaning unverified (GK-E) |  | as row 12; [X] raw :12212. Meaning: OMX, as row 7. Import: pcbnew/pcb_io/eagle/pcb_io_eagle.cpp:175-224 |
| 15 | `mdViaViaSameLayer` | length | clearance.min [board_floor, type(via) vs type(via), pair same_layer] | lossy: exact scope unclear; meaning unverified (GK-E) |  | golden-rules/eez-psu-h24005.txt:48; golden-rules/moteus.txt:91, :115 |
| 16 | `mnLayersViaInSmd` | count | vias_under_smd (provisional) | lossy: meaning unknown (Q-E2) |  | golden-rules/eez-psu-h24005.txt:48; golden-rules/moteus.txt:91, :115 |
| 17 | `mdCopperDimension` | length | edge_clearance.min [board_floor] | exact | meaning from OMX | golden-rules/eez-psu-h24005.txt:48; golden-rules/moteus.txt:91 |
| 18 | `mdDrill` | length | hole_to_hole.min [board_floor] | lossy: meaning from OMX only; see Q-E3 |  | golden-rules/eez-psu-h24005.txt:48; golden-rules/moteus.txt:91 |
| 19 | `mdSmdStop` | length | none | none: meaning unknown (Q-E2) |  | golden-rules/eez-psu-h24005.txt:48; golden-rules/moteus.txt:91 |
| 20 | `msWidth` | length | track_width.min [board_floor] | exact | meaning from OMX | golden-rules/eez-psu-h24005.txt:48; golden-rules/moteus.txt:91 |
| 21 | `msDrill` | length | hole_size.min [board_floor] | exact | meaning from OMX | golden-rules/eez-psu-h24005.txt:48; golden-rules/moteus.txt:91, :115 |
| 22 | `msMicroVia` | length | hole_size.min [board_floor, type(via.micro)] (+ via_span when above every drill) | lossy: disabling reading inferred |  | golden-rules/eez-psu-h24005.txt:48; golden-rules/moteus.txt:91 |
| 23 | `msBlindViaRatio` | ratio | hole_aspect_ratio.max = 1 / ratio [type(via.blind)] | lossy: inverted ratio; meaning unverified (GK-E) |  | golden-rules/eez-psu-h24005.txt:48; golden-rules/moteus.txt:91, :115 |
| 24 | `rvPadTop` | fraction | annular_ratio.ratio [type(pad.pth), layers cu:1] | exact |  | golden-rules/eez-psu-h24005.txt:48; golden-rules/moteus.txt:91, :115 |
| 25 | `rvPadInner` | fraction | annular_ratio.ratio [type(pad.pth), layers inner] | exact |  | as row 24; [X] raw :12223. Import: pcbnew/pcb_io/eagle/pcb_io_eagle.cpp:175-224 |
| 26 | `rvPadBottom` | fraction | annular_ratio.ratio [type(pad.pth), layers cu:N] | exact |  | as row 24; [X] raw :12224. Import: pcbnew/pcb_io/eagle/pcb_io_eagle.h:71, :116 |
| 27 | `rvViaOuter` | fraction | annular_ratio.ratio [type(via), layers outer] | exact |  | as row 24; [X] raw :12225. Meaning: pcbnew/pcb_io/eagle/pcb_io_eagle.h:125. Import: pcbnew/pcb_io/eagle/pcb_io_eagle.cpp:217-218, :2890-2916 |
| 28 | `rvViaInner` | fraction | annular_ratio.ratio [type(via), layers inner] | exact |  | as row 24; [X] raw :12226. Import: pcbnew/pcb_io/eagle/pcb_io_eagle.cpp:2893, :2900 |
| 29 | `rvMicroViaOuter` | fraction | annular_ratio.ratio [type(via.micro), layers outer] | exact |  | as row 24; [X] raw :12227. Import: pcbnew/pcb_io/eagle/pcb_io_eagle.cpp:175-224 |
| 30 | `rvMicroViaInner` | fraction | annular_ratio.ratio [type(via.micro), layers inner] | exact |  | as row 24; [X] raw :12228. Import: pcbnew/pcb_io/eagle/pcb_io_eagle.cpp:175-224 |
| 31 | `rlMinPadTop` | length | annular_width.min + annular_ratio.min [type(pad.pth), layers cu:1] | exact |  | golden-rules/eez-psu-h24005.txt:48; golden-rules/moteus.txt:91 |
| 32 | `rlMaxPadTop` | length | annular_width.max + annular_ratio.max [type(pad.pth), layers cu:1] | exact |  | as row 31; [X] raw :12230. Meaning: OMX: "Rest ring max." |
| 33 | `rlMinPadInner` | length | annular_width.min + annular_ratio.min [type(pad.pth), layers inner] | exact |  | as row 31; [X] raw :12231. Import: pcbnew/pcb_io/eagle/pcb_io_eagle.cpp:209-210, :1862-1864, :1875-1889 |
| 34 | `rlMaxPadInner` | length | annular_width.max + annular_ratio.max [type(pad.pth), layers inner] | exact |  | as row 32; [X] raw :12232. Import: pcbnew/pcb_io/eagle/pcb_io_eagle.cpp:211-212 |
| 35 | `rlMinPadBottom` | length | annular_width.min + annular_ratio.min [type(pad.pth), layers cu:N] | exact |  | as row 31; [X] raw :12233. Import: pcbnew/pcb_io/eagle/pcb_io_eagle.cpp:213-214, :1862-1873 |
| 36 | `rlMaxPadBottom` | length | annular_width.max + annular_ratio.max [type(pad.pth), layers cu:N] | exact |  | as row 32; [X] raw :12234. Import: pcbnew/pcb_io/eagle/pcb_io_eagle.cpp:215-216 |
| 37 | `rlMinViaOuter` | length | annular_width.min + annular_ratio.min [type(via), layers outer] | exact |  | as row 31; [X] raw :12235. Meaning: pcbnew/pcb_io/eagle/pcb_io_eagle.h:126. Import: pcbnew/pcb_io/eagle/pcb_io_eagle.cpp:219-220, :2898, :2912-2914 |
| 38 | `rlMaxViaOuter` | length | annular_width.max + annular_ratio.max [type(via), layers outer] | exact |  | as row 32; [X] raw :12236. Meaning: pcbnew/pcb_io/eagle/pcb_io_eagle.h:127. Import: pcbnew/pcb_io/eagle/pcb_io_eagle.cpp:221-222, :2898, :2912-2914 |
| 39 | `rlMinViaInner` | length | annular_width.min + annular_ratio.min [type(via), layers inner] | exact |  | as row 31; [X] raw :12237. Import: pcbnew/pcb_io/eagle/pcb_io_eagle.cpp:175-224 |
| 40 | `rlMaxViaInner` | length | annular_width.max + annular_ratio.max [type(via), layers inner] | exact |  | as row 32; [X] raw :12238. Import: pcbnew/pcb_io/eagle/pcb_io_eagle.cpp:175-224 |
| 41 | `rlMinMicroViaOuter` | length | annular_width.min + annular_ratio.min [type(via.micro), layers outer] | exact |  | golden-rules/eez-psu-h24005.txt:48; golden-rules/moteus.txt:91 |
| 42 | `rlMaxMicroViaOuter` | length | annular_width.max + annular_ratio.max [type(via.micro), layers outer] | exact |  | as row 41; [X] raw :12240. Import: pcbnew/pcb_io/eagle/pcb_io_eagle.cpp:175-224 |
| 43 | `rlMinMicroViaInner` | length | annular_width.min + annular_ratio.min [type(via.micro), layers inner] | exact |  | as row 41; [X] raw :12241. Import: pcbnew/pcb_io/eagle/pcb_io_eagle.cpp:175-224 |
| 44 | `rlMaxMicroViaInner` | length | annular_width.max + annular_ratio.max [type(via.micro), layers inner] | exact |  | as row 41; [X] raw :12242. Import: pcbnew/pcb_io/eagle/pcb_io_eagle.cpp:175-224 |
| 45 | `useDiameter` | integer bit mask | annular_ratio.honour_library_diameter | lossy: bit layout unknown (Q-E2) |  | golden-rules/eez-psu-h24005.txt:48; golden-rules/moteus.txt:91, :115 |
| 46 | `psTop` | integer enum | pad_shape_rule.shape_top | lossy: shape codes unconfirmed (Q-E4) |  | golden-rules/eez-psu-h24005.txt:48; golden-rules/moteus.txt:91 |
| 47 | `psBottom` | integer enum | pad_shape_rule.shape_bottom | lossy: shape codes unconfirmed (Q-E4) |  | as row 46; [X] raw :12244. Meaning: pcbnew/pcb_io/eagle/pcb_io_eagle.h:104. Import: pcbnew/pcb_io/eagle/pcb_io_eagle.cpp:199-200, :2031-2032 |
| 48 | `psFirst` | integer enum | pad_shape_rule.shape_first | lossy: shape codes unconfirmed (Q-E4) |  | golden-rules/eez-psu-h24005.txt:48; golden-rules/moteus.txt:91, :115 |
| 49 | `psElongationLong` | integer % | pad_shape_rule.elongation_long (%) | exact |  | as row 46; [X] raw :12246. Meaning: pcbnew/pcb_io/eagle/pcb_io_eagle.h:87-89 |
| 50 | `psElongationOffset` | integer % | pad_shape_rule.elongation_offset (%) | exact |  | as row 46; [X] raw :12247. Meaning: pcbnew/pcb_io/eagle/pcb_io_eagle.h:91 |
| 51 | `srRoundness` | fraction | pad_shape_rule.roundness (ratio) | exact |  | as row 46; [X] raw :12255. Meaning: pcbnew/pcb_io/eagle/pcb_io_eagle.h:107 |
| 52 | `srMinRoundness` | length | pad_shape_rule.roundness_min | exact |  | as row 46; [X] raw :12256. Meaning: pcbnew/pcb_io/eagle/pcb_io_eagle.h:109-110. Import: pcbnew/pcb_io/eagle/pcb_io_eagle.cpp:193-194, :2563 |
| 53 | `srMaxRoundness` | length | pad_shape_rule.roundness_max | exact |  | as row 46; [X] raw :12257. Meaning: pcbnew/pcb_io/eagle/pcb_io_eagle.h:112-113. Import: pcbnew/pcb_io/eagle/pcb_io_eagle.cpp:195-196, :2565 |
| 54 | `slThermalIsolate` | length | thermal_relief_gap.min | exact | meaning from OMX | golden-rules/eez-psu-h24005.txt:48; golden-rules/moteus.txt:91 |
| 55 | `slThermalsForVias` | flag | zone_connection.value [type(via)] (0 -> solid) | lossy: meaning unverified (GK-E) |  | as row 54; [X] raw :12259. Meaning: OMX: "Whether to create thermals. 0/1." Import: pcbnew/pcb_io/eagle/pcb_io_eagle.cpp:2871-2945 |
| 56 | `mvStopFrame` | fraction | solder_mask_expansion.ratio [pads and vias] | exact |  | golden-rules/eez-psu-h24005.txt:48; golden-rules/moteus.txt:91 |
| 57 | `mvCreamFrame` | fraction | solder_paste_margin.ratio (negative) [type(pad.smd)] | exact |  | as row 56; [X] raw :12249. Meaning: pcbnew/pcb_io/eagle/pcb_io_eagle.h:96-97 |
| 58 | `mlMinStopFrame` | length | solder_mask_expansion.min clamp | exact |  | as row 56; [X] raw :12250. Meaning: pcbnew/pcb_io/eagle/pcb_io_eagle.h:98. Import: pcbnew/pcb_io/eagle/pcb_io_eagle.cpp:183-184, :2618 |
| 59 | `mlMaxStopFrame` | length | solder_mask_expansion.max clamp | exact |  | as row 56; [X] raw :12251. Meaning: pcbnew/pcb_io/eagle/pcb_io_eagle.h:99. Import: pcbnew/pcb_io/eagle/pcb_io_eagle.cpp:185-186, :2620 |
| 60 | `mlMinCreamFrame` | length | solder_paste_margin.clamp_min | exact |  | as row 56; [X] raw :12252. Meaning: pcbnew/pcb_io/eagle/pcb_io_eagle.h:100. Import: pcbnew/pcb_io/eagle/pcb_io_eagle.cpp:187-188, :2585 |
| 61 | `mlMaxCreamFrame` | length | solder_paste_margin.clamp_max | exact |  | as row 56; [X] raw :12253. Meaning: pcbnew/pcb_io/eagle/pcb_io_eagle.h:101. Import: pcbnew/pcb_io/eagle/pcb_io_eagle.cpp:189-190, :2587 |
| 62 | `mlViaStopLimit` | length | via_protection.tent_max_drill | exact | meaning from OMX | golden-rules/eez-psu-h24005.txt:48; golden-rules/moteus.txt:91 |
| 63 | `dpMaxLengthDifference` | length | skew.max with within_diff_pairs=true | lossy: meaning unverified (GK-E) |  | golden-rules/eez-psu-h24005.txt:48; golden-rules/moteus.txt:91 |
| 64 | `dpGapFactor` | ratio | diff_pair_gap.opt = factor x base | lossy: base unknown (Q-E2) |  | golden-rules/eez-psu-h24005.txt:48; golden-rules/moteus.txt:91, :115 |
| 65 | `checkGrid` | flag | board_off_grid @enabled | lossy: meaning unverified (GK-E) |  | golden-rules/eez-psu-h24005.txt:48; golden-rules/moteus.txt:103 |
| 66 | `checkAngle` | flag | routing_corners (style 45) @enabled | lossy: meaning unverified (GK-E) |  | golden-rules/eez-psu-h24005.txt:48; golden-rules/moteus.txt:91 |
| 67 | `checkFont` | flag | text_font @enabled | lossy: meaning unverified (GK-E) |  | as row 66; [X] raw :12264. Meaning: unverified (general knowledge |
| 68 | `checkRestrict` | flag | disallow @enabled for restrict-layer areas | lossy: meaning unverified (GK-E) |  | golden-rules/eez-psu-h24005.txt:48; golden-rules/moteus.txt:91 |
| 69 | `checkStop` | flag | solder_mask_bridge @enabled | lossy: meaning unverified (GK-E) |  | golden-rules/moteus.txt:91, :107, :115. Meaning: unverified (GK-E). Import: pcbnew/pcb_io/eagle/pcb_io_eagle.cpp:175-224 |
| 70 | `checkValues` | flag | part_value_present @enabled | lossy: meaning unverified (GK-E) |  | golden-rules/moteus.txt:91, :115. Meaning: unverified (GK-E). Import: pcbnew/pcb_io/eagle/pcb_io_eagle.cpp:175-224 |
| 71 | `checkNames` | flag | erc_annotation[missing_name, duplicate_reference] @enabled | lossy: meaning unverified (GK-E) |  | golden-rules/moteus.txt:91, :115. Meaning: unverified (GK-E). Import: pcbnew/pcb_io/eagle/pcb_io_eagle.cpp:175-224 |
| 72 | `checkWireStubs` | flag | dangling_track @enabled | lossy: meaning unverified (GK-E) |  | golden-rules/moteus.txt:91, :115. Meaning: unverified (GK-E). Import: pcbnew/pcb_io/eagle/pcb_io_eagle.cpp:175-224 |
| 73 | `checkPolygonWidth` | flag | zone_min_width @enabled (check) | lossy: meaning unverified (GK-E) |  | golden-rules/moteus.txt:91, :115. Meaning: unverified (GK-E). Import: pcbnew/pcb_io/eagle/pcb_io_eagle.cpp:175-224 |
| 74 | `maxErrors` | count | drc_engine_settings.max_violations | exact |  | golden-rules/eez-psu-h24005.txt:48; golden-rules/moteus.txt:91, :115 |

#### 4.3.3 Net classes and the clearance matrix

| # | Field | Canonical target | Quality | Note | Provenance |
|---|---|---|---|---|---|
| 1 | `class/@number` | @native class number (key of the class map) | model |  | common/io/eagle/eagle_parser.h:1492; common/io/eagle/eagle_parser.cpp:1882 |
| 2 | `class/@name` | @netclass name; 'default' -> Default | model |  | common/io/eagle/eagle_parser.cpp:1883; pcbnew/pcb_io/eagle/pcb_io_eagle.cpp:2665-2673 |
| 3 | `class/@width` | track_width.min [netclass] | exact | a minimum, not KiCad's opt | common/io/eagle/eagle_parser.h:1494; common/io/eagle/eagle_parser.cpp:1884 |
| 4 | `class/@drill` | hole_size.min [netclass, type(via)] | exact | a minimum | common/io/eagle/eagle_parser.h:1495; common/io/eagle/eagle_parser.cpp:1885 |
| 5 | (no Eagle field) via diameter | none | none: Eagle classes have no via-diameter field |  | pcbnew/pcb_io/eagle/pcb_io_eagle.cpp:2676, :2918-2922, :434-435; common/io/eagle/eagle_parser.h:1489-1497 |
| 6 | `<clearance class="0">` inside class N, N ≠ 0 | clearance.min [netclass_exact(N) vs netclass_exact(0)] | exact |  | pcbnew/pcb_io/eagle/pcb_io_eagle.cpp:2682-2688; pcbnew/drc/drc_engine.cpp:456-468 |
| 7 | `<clearance class="0">` inside class 0 | clearance.min [netclass_exact(0) vs netclass_exact(0)] | exact |  | pcbnew/pcb_io/eagle/pcb_io_eagle.cpp:2682-2688; values: golden-rules/eez-psu-h24005.txt:40, :44 |
| 8 | `<clearance class="M">` inside class N, M ≠ 0 (M = N is the class's clearance to itself) | clearance.min [netclass_exact(N) vs netclass_exact(M)] | exact |  | pcbnew/pcb_io/eagle/pcb_io_eagle.cpp:2694-2719; values: golden-rules/eez-psu-h24005.txt:41-42, :50 |
| 9 | No clearance entry for a pair | @netclass: no pair rule; the md* rules apply | model | Eagle's fallback unverified (Q-E5) | pcbnew/pcb_io/eagle/pcb_io_eagle.cpp:421-429; include/netclass.h:141 |
| 10 | Rules sidecar | none | none: KiCad importer artefact (sidecar file) |  | pcbnew/pcb_io/eagle/pcb_io_eagle.cpp:2694, :451-454; pcbnew/pcb_io/eagle/pcb_io_eagle.h:155 |
| 11 | `signal/@class` | @netclass membership | model |  | pcbnew/pcb_io/eagle/pcb_io_eagle.cpp:2779-2789; common/io/eagle/eagle_bin_parser.cpp:3069-3070 |
| 12 | Stale-rules edge case | none | none: KiCad importer defect, not Eagle data |  | pcbnew/pcb_io/eagle/pcb_io_eagle.cpp:2648-2649, :2694, :493-509, :453-454 |

Golden classes, canonical import (rows 3, 4, 7, 8 above): EEZ consolidated and panel `default`, `dgnd`, `gnd` each
get `track_width.min` 0.254 mm and `hole_size.min` 0.5842 mm per class, and clearance rules default-default,
dgnd-dgnd and gnd-gnd of 0.254 mm (golden-rules/eez-psu-h24005.txt:40-46); moteus `default` gets
`track_width.min` 0.127 mm and `CAN` 0.2159 mm, plus a CAN-CAN clearance rule of 0.127 mm (golden-rules/moteus.txt:92-93),
which P8 raises to the 0.1524 mm `mdWireWire` floor unless Q-E5 shows Eagle does otherwise. KiCad's importer instead
gives every class without a class-0 entry the 0.2 mm netclass default, drops the class widths and drills, and leaves the
CAN-CAN rule below the board minimum (inferred from pcbnew/pcb_io/eagle/pcb_io_eagle.cpp:421-432, :2694-2719 and
pcbnew/drc/drc_engine.cpp:2025-2026).

#### 4.3.4 Per-object attributes that act as rules

| # | Object attribute | Canonical target | Quality | Note | Provenance |
|---|---|---|---|---|---|
| 1 | polygon `isolate` | clearance.min [item_override on the zone] | exact |  | common/io/eagle/eagle_parser.h:972, :978; pcbnew/pcb_io/eagle/pcb_io_eagle.cpp:1697-1700 |
| 2 | polygon `width` | zone_min_width.min (pen width) | lossy: the pen also inflates the outline |  | common/io/eagle/eagle_parser.h:968; pcbnew/pcb_io/eagle/pcb_io_eagle.cpp:1654-1657, :1692-1695 |
| 3 | polygon `spacing` | zone_fill_style.hatch_gap = spacing - width | exact |  | common/io/eagle/eagle_parser.h:970; pcbnew/pcb_io/eagle/pcb_io_eagle.cpp:1682-1690 |
| 4 | polygon `pour` | zone_fill_style.fill_mode; cutout -> disallow zone_fill | exact |  | common/io/eagle/eagle_parser.h:971; pcbnew/pcb_io/eagle/pcb_io_eagle.cpp:1671-1681 |
| 5 | polygon `thermals` | zone_connection [zone] | exact |  | common/io/eagle/eagle_parser.h:974; pcbnew/pcb_io/eagle/pcb_io_eagle.cpp:1703-1714 |
| 6 | polygon `orphans` | zone_fill_style.island_removal | exact |  | common/io/eagle/eagle_parser.h:973; pcbnew/pcb_io/eagle/pcb_io_eagle.cpp:1605-1610 |
| 7 | polygon `rank` | zone_fill_style.priority = 6 - rank | exact |  | common/io/eagle/eagle_parser.h:975, :981, :987-992; pcbnew/pcb_io/eagle/pcb_io_eagle.cpp:1716-1717 |
| 8 | pad or smd `thermals="no"` | zone_connection = solid [item_override] | exact |  | common/io/eagle/eagle_parser.h:855, :914; pcbnew/pcb_io/eagle/pcb_io_eagle.cpp:2622-2624 |
| 9 | pad `stop`; smd `stop` and `cream` | none (object data: pad mask and paste layers) | none: object data, not a rule |  | common/io/eagle/eagle_parser.h:854, :913, :915; pcbnew/pcb_io/eagle/pcb_io_eagle.cpp:2040-2042, :2589-2605 |
| 10 | smd `roundness` | pad_shape_rule.roundness [item] | exact |  | common/io/eagle/eagle_parser.h:911, :921; pcbnew/pcb_io/eagle/pcb_io_eagle.cpp:2571-2576 |
| 11 | pad `diameter`, via `diameter` | none (object geometry) | none: object data; annular rules still apply |  | common/io/eagle/eagle_parser.h:560, :880; pcbnew/pcb_io/eagle/pcb_io_eagle.cpp:2079-2083, :2890-2916, :1855-1860 |
| 12 | via `alwaysstop` | via_protection (opening forced) [item_override] | exact |  | common/io/eagle/eagle_parser.h:562, :573; common/io/eagle/eagle_parser.cpp:779-784 |
| 13 | tRestrict, bRestrict (41, 42); vRestrict (43) | disallow [in_area(shape)]: t/bRestrict -> track, via, zone_fill, pad on cu:1 / cu:N; vRestrict -> via on all copper | exact |  | common/io/eagle/eagle_parser.h:1240-1242; pcbnew/pcb_io/eagle/pcb_io_eagle.cpp:121-147 |

#### 4.3.5 Eagle-specific semantics and how the canonical model handles them

| # | Eagle semantics | Canonical handling | Provenance |
|---|---|---|---|
| 1 | Restring percentage with min and max | annular_ratio with min and max clamps | pcbnew/pcb_io/eagle/pcb_io_eagle.cpp:2084-2091, :2895-2900; values: golden-rules/moteus.txt:91 |
| 2 | Min greater than max | clamp rule: result = min(max, max(min, value)); canonical keeps both bounds and applies the same order | pcbnew/pcb_io/eagle/pcb_io_eagle.cpp:102-109 |
| 3 | Restring per layer | one annular rule per layer group with a layer filter | pcbnew/pcb_io/eagle/pcb_io_eagle.cpp:2082, :2090, :2893, :2900, :1862-1864, :1875-1889 |
| 4 | Library diameter against restring | annular_ratio.honour_library_diameter (bits unresolved) | pcbnew/pcb_io/eagle/pcb_io_eagle.h:301-317; pcbnew/pcb_io/eagle/pcb_io_eagle.cpp:2079-2083, :1840-1891 |
| 5 | Stop-mask percentage | solder_mask_expansion.ratio with clamps | pcbnew/pcb_io/eagle/pcb_io_eagle.cpp:2615-2620; pcbnew/drc/drc_engine.cpp:1312-1340 |
| 6 | Cream percentage | solder_paste_margin.ratio (negative) with clamps | pcbnew/pcb_io/eagle/pcb_io_eagle.cpp:2585-2587 |
| 7 | Via stop limit | via_protection.tent_max_drill; via alwaysstop as item override | pcbnew/board_design_settings.cpp:242-243; OMX, `mlViaStopLimit` entry |
| 8 | Layer-setup string syntax | stackup + via_span (syntax unverified) | golden-rules/moteus.txt:115; [P] raw :14534 |
| 9 | `mtCopper` and `mtIsolate` indexing | stackup thickness lists (indexing hypothesis) | golden-rules/moteus.txt:91; golden-rules/eez-psu-h24005.txt:48 |
| 10 | Unit suffixes | units: any value without 'mil' is mm; canonical stores mm | pcbnew/pcb_io/eagle/pcb_io_eagle.cpp:88-99; common/io/eagle/eagle_parser.cpp:219-261 |
| 11 | Version-dependent parameter sets | absent parameters are absent rules (no default injected) | golden-rules/eez-psu-h24005.txt:48; golden-rules/moteus.txt:91, :103 |
| 12 | ERULES defaults (missing block or parameter, and every `.lbr` load) | tool defaults are not imported as rules (origin tool_default only if written) | pcbnew/pcb_io/eagle/pcb_io_eagle.h:51-83, :312-313; pcbnew/pcb_io/eagle/pcb_io_eagle.cpp:493-509, :3330-3345 |
| 13 | Eagle 4/5 binary boards | binary boards: 5 recoverable values only | common/io/eagle/eagle_bin_parser.cpp:10, :1286-1346, :2107-2121; common/io/eagle/eagle_bin_parser.h:117-124 |
| 14 | Pairwise matrix against per-class clearance | binary clearance rules per class pair (symmetric) | pcbnew/drc/drc_engine.cpp:456-468; pcbnew/drc/drc_rule_condition.cpp:105-117 |
| 15 | Explicit pair rules are not floored (inferred) | floors: canonical floors apply unless floor_exempt (section 5, P8) | pcbnew/drc/drc_engine.cpp:2025-2026, :2089-2160 |
| 16 | Zero values | zero values kept as written; 0 mil md* flagged as 'unchecked?' (Q-E5) | golden-rules/moteus.txt:91; common/io/eagle/eagle_parser.h:1494-1495 |

### 4.4 Export recipes (canonical -> tool)

**To KiCad.** Write one `.kicad_dru` file plus board-setup values:
1. `board_floor` rules of `clearance`, `track_width`, `connection_width`, `via_diameter`, `hole_size`,
   `annular_width`, `hole_to_hole`, `hole_clearance`, `edge_clearance`, `silk_clearance`, `text_height`,
   `text_thickness`, `starved_thermal` go to `board.design_settings.rules.min_*` (section 4.1.3). Where a canonical
   floor has a non-`all` scope, write it as an explicit DRU rule instead.
2. `netclass` rules go to the netclass fields (section 4.1.6) when they have a netclass-only scope and set the field
   KiCad stores (clearance min; width, via and drill opt; pair width and gap opt). Everything else becomes an explicit
   rule.
3. Explicit rules become one DRU rule each, written from the largest priority number down to 1, so the priority-1
   rule comes last and wins under KiCad's last-match resolution (section 5, P3). Rules that are not `floor_exempt` and
   have a min below a floor are written with min = the floor, because KiCad would not floor them (P8). Shape:

   ```
   (rule "<native.name or id>"
     (constraint <keyword> (min <v><unit>) (opt <v><unit>) (max <v><unit>))
     (condition "<scope expression translated back to KiCad syntax>")
     (layer <layer>)
     (severity <severity>))
   ```
4. Built-in checks with no constraint (DRC codes marked "lossy: board-wide severity only") export only their
   severity, and only when the canonical rule is unscoped. A scoped `short_circuit.allowed=true` (for example an
   Altium ShortCircuit rule with `ALLOWED=TRUE` on a net-tie scope) is approximated by a clearance rule with
   `(min 0mm)` on that scope, because KiCad tests shorts and track crossings only when the pair's clearance is > 0 and
   not ignored (pcbnew/drc/drc_test_provider_copper_clearance.cpp:275-281, :866-879); the approximation also removes the
   clearance check for that scope, so the rule is marked lossy.
5. Check the target version (section 4.1.10) and drop or comment out what it cannot read, as KiCad's own downgrade
   filter does (include/drc_rules_downgrade.h:472-535).

**To Altium.** One rule record per canonical rule: `RULEKIND` from section 4.2.1, `NAME` = `native.name` or id,
`PRIORITY` = priority renumbered from 1 within each kind, `SCOPE1EXPRESSION` / `SCOPE2EXPRESSION` from `expr` /
`expr_b` (section 2.5, reverse), `NETSCOPE` / `LAYERKIND` from `scope.pair`, parameters in mil with the `mil` suffix
(KiCad's reader rejects values without it, common/io/altium/altium_props_utils.cpp:86-113). Floors become the lowest
priority `All` / `All` rule of the kind, and every rule of the kind that is not `floor_exempt` is written with
min = max(min, floor), as in the KiCad recipe (step 3): Altium applies only the highest-priority matching rule (section
4.2.3), so a lowest-priority `All` / `All` rule is a default, not a floor, and would never raise a scoped rule's smaller
min. Per-layer canonical rules of one kind fold into the per-layer keys of a
single Width or DiffPairsRouting record. Kinds marked "none" for Altium are not written.

**To Eagle.** One `<designrules>` block: `md*`, `ms*`, `rl*`, `rv*`, `mv*`, `ml*`, `sl*`, `dp*`, `ps*`, `sr*` from the
board-wide canonical rules of the matching kinds; class `width` / `drill` from netclass-scoped minima; the clearance
matrix from class-pair clearance rules. Any canonical rule with a scope Eagle cannot express (component, area, layer
other than the restring groups) is dropped and reported.

## 5. Precedence and conflict resolution

### 5.1 Canonical resolution (P1-P13)

- **P1 Tiers.** Every rule has an `origin`: `fab_floor` (W7 capability), `board_floor` (board-wide minima: KiCad
  `rules.min_*` and their KiCad 5 equivalents, Eagle `ms*` and `md*`; Altium has no separate floor),
  `tool_default` (hard-coded tool defaults written into a file), `netclass`, `item_default` (object-level settings
  that apply only when no explicit rule matches: KiCad's footprint zone connection and a zone's own connection, thermal
  gap and spoke width, P9), `explicit` (DRU rules, scoped Altium rules, Eagle class matrix entries, W8 interface rules)
  and `item_override` (pad and footprint local values that KiCad checks before any rule, and zone local clearance, P9).
- **P2 Candidates.** For a query (kind, object or pair, layer) the candidates are the enabled rules of that kind whose
  layer filter contains the layer under test and whose scope matches. Source order is kept. (KiCad: implicit rules
  first, then file order, pcbnew/drc/drc_engine.cpp:906-927, :744-767.)
- **P3 Selection.** Candidates are ranked by tier, `item_override` > `explicit` > `item_default` > `netclass` >
  `tool_default` > `board_floor` and `fab_floor`, then by the lowest `priority` number, then by the later definition.
  - Numbers: for each of `min`, `opt` and `max` separately, the winner is the top-ranked candidate that sets the
    field. This reproduces KiCad, where a matching rule overwrites min, opt and max only when it sets them and the last
    one wins per value (pcbnew/drc/drc_engine.cpp:1083-1105, :1997-2003).
  - Everything else comes from the top-ranked matching candidate as a whole, whether or not it sets a number: `value`
    tokens (zone connection, disallow kinds), `options` (within_diff_pairs; time or space domain) and the severity
    (P13). This also reproduces KiCad, which copies the disallow flags, the zone connection, the parent rule (and so
    the severity) and the option bits from every later matching rule (pcbnew/drc/drc_engine.cpp:1129-1136;
    `SetOptionsFromOther`, pcbnew/drc/drc_rule.h:198). Example: a later `skew` rule that sets only `(min ...)` and has
    no `within_diff_pairs` clears the option, and its severity applies to the max inherited from an earlier rule.
    When the min, opt and max winners use different unit domains (ps against mm), the kernel reports the mix as a
    finding instead of comparing mixed units (KiCad's parser rejects mixed units only within one rule,
    pcbnew/drc/drc_rule_parser.cpp:491-495).
  - Floors are ordinary candidates at the bottom of the order, so a kind whose only candidates are floors is still
    checked against them, as KiCad checks against its board-setup implicit rules (pcbnew/drc/drc_engine.cpp:163-236).
    With no DRU rules this is the case in KiCad for hole_to_hole, hole_clearance, edge_clearance, connection_width,
    silk_clearance, text_height, text_thickness and starved_thermal. P8 then raises any selected min.
  - Only when no candidate at all matches is there no constraint, and the check is skipped (KiCad `NULL_CONSTRAINT`,
    pcbnew/drc/drc_engine.cpp:2234-2238).
  KiCad's help says matching stops at the first match (pcbnew/dialogs/panel_setup_rules_help_6notes.md:7-9), which the
  code does not do (Q-K3). The model also reproduces Altium, whose winning rule sets every parameter, so per-field
  selection equals Altium's whole-rule first match (AD13).
- **P4 Matching.** Binary scopes are symmetric by default (KiCad tries (A, B) then (B, A),
  pcbnew/drc/drc_rule_condition.cpp:105-117; Altium's behaviour is assumed the same, unverified (GK)). A rule whose
  layer filter excludes the layer under test is skipped (pcbnew/drc/drc_engine.cpp:1819-1835). Exception reproduced
  from KiCad: edge_clearance and the silk-to-edge check are evaluated with no layer under test (`UNDEFINED_LAYER`,
  pcbnew/drc/drc_test_provider_edge_clearance.cpp:108-110, :176, :218), and the engine tests a rule's `(layer ...)` only
  when the query has a board layer; otherwise it only requires the rule's layers to overlap the board's enabled layers
  (pcbnew/drc/drc_engine.cpp:1819-1820). A KiCad rule `(layer F.Cu)` on edge_clearance therefore also governs B.Cu
  copper. Such rules are imported without `layers` (the KiCad layer clause stays in `native.text` and is written back
  on export), so canonical results match KiCad.
- **P5 Netclass resolution.** A net's effective netclass is composed field by field from all classes it belongs to,
  highest-priority class first (lower number wins; Default last), common/project/net_settings.cpp:1141-1165,
  :1192-1240.
- **P6 Two netclasses on one pair.** When A and B belong to different netclasses, the larger of the two netclass
  values wins (KiCad sorts netclass clearance rules ascending so the larger applies last,
  pcbnew/drc/drc_engine.cpp:456-468; fast path max(A, B), :1929-1994).
- **P7 Differential pairs.** The resolved `diff_pair_gap.min` of a pair is raised to the resolved clearance between
  the two tracks (pcbnew/drc/drc_engine.cpp:2005-2022); a netclass whose pair gap is below its clearance gets an
  extra clearance rule for coupled pairs (:353-371).
- **P8 Floors.** After selection, a lower bound is raised to the floors,
  `min := max(min, board_floor.min, fab_floor.min)`, for every kind whose min is a manufacturability or safety bound
  (clearances, widths, sizes, ring, sliver, spacing, text). Floors never lower a value and do not touch `opt` or `max`; an `opt` below a floor is
  reported. A rule with `floor_exempt: true` skips the board floor (never the fab floor without a waiver). Imported
  KiCad DRU rules get `floor_exempt: true`, because KiCad returns an explicit rule's result before applying
  `rules.min_clearance` (pcbnew/drc/drc_engine.cpp:2025-2026 against :2143-2159; inferred, Q-K10). The kernel then
  reports every explicit rule below a floor as a finding; the moteus Eagle class rule CAN-to-CAN 0.127 mm against a
  0.1524 mm board minimum (golden-rules/moteus.txt:91, :93) is the golden example.
- **P9 Item overrides and object defaults.** Per kind, from first to last; the first source that gives a value wins.
  "Rules" means explicit rules: a match by an implicit rule does not stop the search (pcbnew/drc/drc_engine.cpp:2025-2026,
  :2058-2059).
  - Clearance and hole clearance: pad override, else the pad's parent-footprint override [item_override], floored at
    the board minimum (pcbnew/drc/drc_engine.cpp:1180-1247); then rules. Only PAD implements the override (the pad
    accessor returns the pad's padstack clearance, else the footprint's, pcbnew/pad.cpp:1913-1921†; other connected
    items return none, pcbnew/board_connected_item.h:169-172†). A zone's local clearance is not an override: it is
    max-ed with the rule result when no explicit rule matched (pcbnew/drc/drc_engine.cpp:2089-2160).
  - Zone connection: pad override [item_override] > rules on the pad > rules on the parent footprint > the footprint's
    own zone connection [item_default] > the zone's pad-connection setting [item_default] (pcbnew/drc/drc_engine.cpp:1251-1266,
    :2031-2060, :2164-2197). The footprint setting therefore ranks below DRU rules, while the pad setting ranks above them.
  - Thermal relief gap and spoke width: pad override [item_override] (spoke width raised to the zone's minimum width) >
    rules on the pad > rules on the parent footprint > the zone's own value [item_default] (:1268-1300, :2031-2060,
    :2199-2230).
  - Solder mask expansion, paste margin and paste ratio: pad override (for mask also a shape or track override) > the
    parent footprint's override [both item_override, checked before rules] > rules on the item > rules on the parent
    footprint > board setup [board default] (:1312-1399, :2031-2083). The board-setup fallback is reached only for
    pads in footprints (:2031).
- **P10 Net ties.** Clearance between items joined by a footprint net tie is 0 mm before any rule
  (pcbnew/drc/drc_engine.cpp:1142-1176).
- **P11 Caps.** Clearance-family minima are capped at KiCad's `MAXIMUM_CLEARANCE` on export
  (pcbnew/drc/drc_engine.cpp:1110-1120; value not in the checkout).
- **P12 All-match kinds.** `custom_assertion` tests every matching rule (pcbnew/drc/drc_engine.cpp:2278-2347);
  `disallow` fires on any matching rule (:1716-1817); `room_confinement` obeys every matching rule, because Altium
  states that a component scoped by several Room Definition rules obeys all of them (AD15; verify in L1), and KiCad's
  form of it is an `assertion` (section 4.1.1 row 41), which is all-match already; `length`, `skew`, `via_count` group
  objects by the rule that won for them, so skew compares only nets matched by the same rule
  (pcbnew/drc/drc_test_provider_matched_length.cpp:393-405).
- **P13 Severity.** The severity of the top-ranked matching candidate of P3, taken as a whole, even when the min, opt or
  max that fails came from another candidate; `inherit` resolves to the rule set's default for the kind (section 2.7).
  `ignore` suppresses findings for the matched objects; the tag `waive` emits them as dismissed.

### 5.2 Tool behaviours the importers must preserve

| Topic | KiCad | Altium | Eagle |
|---|---|---|---|
| Order | implicit rules, then DRU file order; last match per value for min, opt and max, last match as a whole for options, tokens and severity (P3) | `PRIORITY` per kind, 1 highest, first match, independent of scope specificity (AD13), except Room Definition, where every matching rule applies (AD15, P12; verify in L1); legacy `TScopeKind` precedence from Board (lowest) to Region (highest), use alongside `PRIORITY` unverified (Q-A7) | one rule set, no order |
| Netclass vs rule | netclass values are implicit rules, so any matching DRU rule overrides them (pcbnew/drc/drc_engine.cpp:262-454, :2025-2026) | a netclass is just a scope (`InNetClass`) of an ordinary prioritised rule | class `width` and `drill` are minima; class matrix entries are pair rules |
| Board minimums | `rules.min_*` are implicit rules and the only candidates when no other rule matches (P3); `min_clearance` and `min_hole_clearance` also floor local overrides and implicit results, not explicit rules (P8); `min_via_diameter` below drill + 2 x annulus has no effect (validator, pcbnew/board_design_settings.cpp:1760-1773) | no separate floor; the lowest-priority `All`/`All` rule is a default, not a floor: a higher-priority scoped rule with a smaller min wins outright (section 4.4 writes floors into every non-exempt rule) | `md*`, `ms*` act as board-wide minima; how they combine with the class matrix is unverified (Q-E5) |
| Layer filter | `(layer ...)` narrows every check except edge_clearance and silk-to-edge, which KiCad evaluates with no layer under test; for those the filter only has to overlap the enabled layers (P4) | per-layer parameter keys and `OnLayer` scopes (section 2.5) | layer groups of the restring and annulus rows (section 2.5) |
| Local values | pad and footprint overrides before rules; footprint zone connection and zone thermal and connection settings after rules; zone local clearance max-ed with the rule result (P9) | pad and via mask modes None / Rule / Manual (pcbnew/pcb_io/altium/altium_parser_pcb.h:183-189) | per-object attributes (section 4.3.4) |
| Disabling | `(severity ignore)` rules; severity map; a clearance rule with min <= 0 mm or severity ignore also turns off short and crossing checks for its scope (section 3.1) | `ENABLED`, DRC check lists, negative `GAP` (ODrive ComponentClearance -100 mm for net ties, golden-rules/odrive-v3.txt:23) | `check*` switches |

Canonical import defaults where the tools are silent or disagree: Eagle `md*` and `ms*` are imported as
`board_floor` and class-matrix entries as `explicit` without `floor_exempt`, so the larger value wins (pending Q-E5);
Altium `All`/`All` rules are imported as `explicit` with the lowest priority of their kind (they stay rules, so an
Altium round trip is exact), and the kernel derives floors from W7 only; a negative Altium clearance (`GAP` < 0) is kept
as written and treated as "check disabled for this scope", the way the ODrive net-tie rules use it.

## 6. Golden-board coverage

### 6.1 Constraint kinds used by the ten golden boards

Every row of the golden-usage union table (80 constraint kinds against boards) maps to the canonical model:
80 rows: 75 map to canonical kinds, 4 to rule-model fields, 1 provisional. Board codes: BMS bms-c1, EEZ eez-psu-h24005, EVB esp32-evb, GLA glasgow, HRF hackrf-one,
JET jetson-nano-baseboard, MOT moteus, ODR odrive-v3, SHE sh-esp32, SPD speeduino-v04.

| # | Constraint kind | Native keys (abridged) | Canonical target | Status | Boards with an active value (abridged) | Provenance |
|---|---|---|---|---|---|---|
| 1 | Copper clearance floor (board minimum) | rules.min_clearance | clearance.min [board_floor] | kind | BMS 0.2 mm, GLA 0.09 mm (revD), HRF 0.127 mm, MOT 0.127 mm | golden-rules/{bms-c1,glasgow,hackrf-one,moteus,esp32-evb,jetson-nano-baseboard,sh-esp32,speeduino-v04}.txt |
| 2 | Net-to-net copper clearance, default or per class | netclass clearance; KiCad 5 trace_clearance, (net_class (clearance)); .pro Clearance; Eagle mdWireWire, mdWirePad, mdWireVia, mdPadPad, mdPadVia, mdViaVia ... | clearance.min [netclass / explicit] | kind | all 10 (range 0.1 mm GLA revD to 1.5 mm EEZ Mains) | golden-rules/odrive-v3.txt:35 |
| 3 | Clearance between two named classes | Eagle &lt;class>&lt;clearance class=...> | clearance [netclass_exact(N) vs netclass_exact(M)] | kind | EEZ (0.254 mm), MOT (CAN–CAN 0.127 mm) | golden-rules/eez-psu-h24005.txt:40-46; golden-rules/moteus.txt:93 |
| 4 | Clearance per object-type pair | Eagle mdSmdPad, mdSmdVia, mdSmdSmd, mdViaViaSameLayer; Altium OBJECTCLEARANCES | clearance [type() vs type()] | kind | EEZ (8 mil), MOT (mdViaViaSameLayer 8 mil), ODR (SMD pad–SMD pad 7.2835 mil) | golden-rules/odrive-v3.txt:35 |
| 5 | Local item clearance override | pad or footprint (clearance) | clearance [item_override] | kind | BMS, EVB, GLA (revA-revC3), HRF, SHE | golden-rules/*.txt (per-board values; keys resolved in section 6.2) |
| 6 | Zone (pour) clearance | zone (connect_pads (clearance)); KiCad 5 zone_clearance; Altium PlaneClearance | clearance [item_override on zones]; plane_clearance (Altium) | kind | all 10 | golden-rules/odrive-v3.txt:33 |
| 7 | Copper-to-hole clearance | rules.min_hole_clearance; Eagle mdDrill | hole_clearance (Eagle `mdDrill` mapped to hole_to_hole, Q-E3) | kind | BMS, EEZ, EVB, GLA, HRF, JET, MOT, SHE, SPD | golden-rules/*.txt (per-board values; keys resolved in section 6.2) |
| 8 | Hole-to-hole spacing | rules.min_hole_to_hole; .pro MinHoleToHole; Altium HoleToHoleClearance | hole_to_hole | kind | BMS, EVB, GLA, HRF, JET, MOT, ODR, SHE, SPD | golden-rules/odrive-v3.txt:14 |
| 9 | Copper to board edge | rules.min_copper_edge_clearance; Eagle mdCopperDimension | edge_clearance | kind | BMS, EEZ, EVB, GLA, JET, MOT, SHE, SPD | golden-rules/*.txt (per-board values; keys resolved in section 6.2) |
| 10 | Minimum track width | rules.min_track_width; KiCad 5 trace_min; .pro MinTrackWidth; Eagle msWidth; Altium Width MINLIMIT | track_width.min [board_floor] | kind | all 10 | golden-rules/odrive-v3.txt:31 |
| 11 | Track width per class or preferred | netclass track_width; Eagle class width; Altium Width PREFEREDWIDTH | track_width.opt [netclass]; Eagle class width -> track_width.min [netclass] | kind | all 10 | golden-rules/*.txt (per-board values; keys resolved in section 6.2) |
| 12 | Maximum track width | Altium Width MAXLIMIT, per-layer MAXWIDTH | track_width.max | kind | ODR | golden-rules/*.txt (per-board values; keys resolved in section 6.2) |
| 13 | Per-layer track width | Altium Width &lt;LAYER>_MINWIDTH/PREFWIDTH/MAXWIDTH | track_width [layers] | kind | ODR | golden-rules/*.txt (per-board values; keys resolved in section 6.2) |
| 14 | Minimum via diameter | rules.min_via_diameter; KiCad 5 via_min_size; .pro MinViaDiameter; Altium RoutingVias MINWIDTH | via_diameter.min [board_floor] | kind | all 10 (via the KiCad files for EEZ and MOT) | golden-rules/odrive-v3.txt:27 |
| 15 | Via size per class or preferred | netclass via_diameter, via_drill; KiCad 5 via_size, via_drill; Eagle class drill; Altium RoutingVias WIDTH, HOLEWIDTH | via_diameter.opt, hole_size.opt [netclass]; Eagle class drill -> hole_size.min | kind | all 10 | golden-rules/odrive-v3.txt:27 |
| 16 | Minimum drill (through hole) | rules.min_through_hole_diameter; KiCad 5 via_min_drill; .pro MinViaDrill; Eagle msDrill; Altium HoleSize MINLIMIT, RoutingVias MINHOLEWIDTH | hole_size.min [board_floor] | kind | all 10 | golden-rules/odrive-v3.txt:20, :27 |
| 17 | Maximum hole size | Altium HoleSize MAXLIMIT; RoutingVias MAXHOLEWIDTH, MAXWIDTH | hole_size.max; via_diameter.max | kind | ODR | golden-rules/odrive-v3.txt:20, :27 |
| 18 | Minimum annular ring | rules.min_via_annular_width; Eagle rlMin* with rv*; Altium MinimumAnnularRing | annular_width.min | kind | all 10 (EEZ through its Eagle files only; its KiCad 5 file has no annulus key) | golden-rules/odrive-v3.txt:3 |
| 19 | Restring maximum and ratio | Eagle rlMax*, rv* | annular_width.max + annular_ratio | kind | EEZ, MOT | golden-rules/*.txt (per-board values; keys resolved in section 6.2) |
| 20 | Microvia minimum size and drill | rules.min_microvia_diameter, min_microvia_drill; KiCad 5 uvia_min_size, uvia_min_drill; .pro MinMicroVia*; Eagle msMicroVia, rlMinMicroVia* | via_diameter.min, hole_size.min [type(via.micro)] | kind | BMS, EEZ, GLA, HRF, JET, MOT, SHE, SPD | golden-rules/*.txt (per-board values; keys resolved in section 6.2) |
| 21 | Microvia size per class | netclass microvia_diameter, microvia_drill; KiCad 5 uvia_size, uvia_drill | via_diameter.opt, hole_size.opt [netclass, type(via.micro)] | kind | all 9 KiCad boards | golden-rules/*.txt (per-board values; keys resolved in section 6.2) |
| 22 | Microvias allowed | rules.allow_microvias (obsolete key); KiCad 5 uvias_allowed; .pro AllowMicroVias; Altium RoutingVias VIASTYLE | via_span.allowed_types | kind | none | golden-rules/odrive-v3.txt:27 |
| 23 | Blind and buried vias | rules.allow_blind_buried_vias (obsolete key); .pro AllowBlindVias; Eagle msBlindViaRatio; Altium LayerPairs | via_span.allowed_types; hole_aspect_ratio (Eagle `msBlindViaRatio`); via_span.enforce_pairs (Altium LayerPairs) | kind | EEZ, MOT (ratio 0.5), ODR (ENFORCE=TRUE) | golden-rules/odrive-v3.txt:21 |
| 24 | Minimum connection width | rules.min_connection | connection_width.min | kind | BMS 0.2 mm | golden-rules/*.txt (per-board values; keys resolved in section 6.2) |
| 25 | Starved thermal (minimum spokes) | rules.min_resolved_spokes; Altium DRC MINSTARVEDCOPPERPERCENT | starved_thermal.min_resolved / .min_copper | kind | BMS 1, EVB 2, GLA 2, JET 2, MOT 2, SHE 2 (count); ODR 50 % | golden-rules/*.txt (per-board values; keys resolved in section 6.2) |
| 26 | Thermal relief geometry | zone thermal_gap, thermal_bridge_width; pad thermal_* overrides; Eagle slThermalIsolate; Altium PlaneConnect, PolygonConnect widths and gaps | thermal_relief_gap + thermal_spoke_width | kind | all 10 | golden-rules/odrive-v3.txt:32, :34 |
| 27 | Zone or plane connection style | zone connect_pads mode; pad zone_connect; Eagle slThermalsForVias; Altium CONNECTSTYLE, PLANECONNECTSTYLE, RELIEFENTRIES | zone_connection (+ thermal_spoke_count for `RELIEFENTRIES`) | kind | all 10 | golden-rules/odrive-v3.txt:32, :34 |
| 28 | Zone minimum width | zone min_thickness | zone_min_width | kind | BMS, EEZ, EVB, GLA, HRF, JET, MOT, SHE, SPD | golden-rules/*.txt (per-board values; keys resolved in section 6.2) |
| 29 | Solder-mask expansion, board | KiCad pad_to_mask_clearance; rules.solder_mask_clearance; .pro SolderMaskClearance; Eagle mvStopFrame, mlMinStopFrame, mlMaxStopFrame; Altium SolderMaskExpansion | solder_mask_expansion [board default] | kind | EEZ, EVB, HRF, MOT (Eagle 2 mil), ODR, SHE (HAT boards), SPD | golden-rules/odrive-v3.txt:25 |
| 30 | Solder-mask expansion, per pad | pad or footprint solder_mask_margin | solder_mask_expansion [item_override] | kind | BMS, EEZ, EVB, GLA, HRF, JET, MOT, SHE | golden-rules/*.txt (per-board values; keys resolved in section 6.2) |
| 31 | Solder-mask minimum web or sliver | KiCad 5 solder_mask_min_width; rules.solder_mask_min_width; .pro SolderMaskMinWidth; Altium MinimumSolderMaskSliver | solder_mask_sliver.min (KiCad values also solder_mask_bridge.min_web, section 4.1.4 row 2) | kind | GLA 0.09 mm, SPD 0.25 mm, ODR 7.0866 mil | golden-rules/odrive-v3.txt:13 |
| 32 | Mask aperture to other-net copper | rules.solder_mask_to_copper_clearance | solder_mask_bridge.min_copper_to_aperture | kind | GLA 0.005 mm (revD) | golden-rules/*.txt (per-board values; keys resolved in section 6.2) |
| 33 | Mask bridges inside footprints | allow_soldermask_bridges_in_footprints | solder_mask_bridge.allow_in_footprint | kind | none | golden-rules/*.txt (per-board values; keys resolved in section 6.2) |
| 34 | Via tenting | KiCad tenting; Eagle mlViaStopLimit | via_protection | kind | BMS, GLA (revD), JET (front and back); EEZ, MOT (25 mil) | golden-rules/*.txt (per-board values; keys resolved in section 6.2) |
| 35 | Paste margin, board | rules.solder_paste_clearance; .pro SolderPasteClearance; Eagle mvCreamFrame, ml*CreamFrame; Altium PasteMaskExpansion | solder_paste_margin.abs [board default] | kind | none | golden-rules/odrive-v3.txt:26 |
| 36 | Paste margin, per pad | pad solder_paste_margin | solder_paste_margin [item_override] | kind | EVB, GLA, JET | golden-rules/*.txt (per-board values; keys resolved in section 6.2) |
| 37 | Paste ratio | pad_to_paste_clearance_ratio; rules.solder_paste_margin_ratio; .pro SolderPasteRatio; pad ratio | solder_paste_margin.ratio | kind | HRF (−0.12 board, pads −0.05 to −0.01) | golden-rules/*.txt (per-board values; keys resolved in section 6.2) |
| 38 | Silk-to-silk clearance | rules.min_silk_clearance; Altium SilkToSilkClearance | silk_clearance | kind | ODR 10 mil | golden-rules/odrive-v3.txt:11 |
| 39 | Silk to mask opening or exposed copper | KiCad severity silk_over_copper with min_silk_clearance; Altium SilkToSolderMaskClearance | silk_to_mask_clearance | kind | ODR 7.874 mil | golden-rules/odrive-v3.txt:12 |
| 40 | Silk to board edge or region | KiCad severity silk_edge_clearance; Altium SilkToBoardRegionClearance | silk_to_edge_clearance | kind | ODR (no value key) | golden-rules/odrive-v3.txt:10 |
| 41 | Silk text minimum height and stroke | rules.min_text_height, min_text_thickness | text_height + text_thickness | kind | BMS, EVB (Rev L), GLA, JET, MOT (0.14 mm / 0.018 mm), SHE | golden-rules/*.txt (per-board values; keys resolved in section 6.2) |
| 42 | Courtyard overlap and component clearance | KiCad severity courtyards_overlap; .pro ProhibitOverlappingCourtyards, RequireCourtyardDefinitions; Altium ComponentClearance | courtyard_clearance | kind | EVB, GLA (revC1), SHE, SPD (.pro 1); GLA (revC3 warning); ODR (3 rules) | golden-rules/odrive-v3.txt:23, :38-39 |
| 43 | Holes inside courtyards | KiCad severities pth_inside_courtyard, npth_inside_courtyard | hole_in_courtyard | kind | JET (warning) | golden-rules/*.txt (per-board values; keys resolved in section 6.2) |
| 44 | Component height | Altium Height | component_height | kind | ODR (0 to 1000 mil, preferred 500 mil) | golden-rules/odrive-v3.txt:22 |
| 45 | Diff-pair width and gap | netclass diff_pair_width, diff_pair_gap; KiCad 5 (net_class (diff_pair_*)); .pro dPair*; Eagle dpGapFactor; Altium DiffPairsRouting | track_width [diffpair] + diff_pair_gap | kind | all 10 (ODR rule has no pairs) | golden-rules/odrive-v3.txt:15 |
| 46 | Diff-pair via gap | netclass diff_pair_via_gap | diff_pair_gap.via_gap | kind | BMS, EVB, GLA, HRF, JET, MOT, SHE, SPD (0.25 mm; router value, no DRC rule per section 4.1.6) | golden-rules/*.txt (per-board values; keys resolved in section 6.2) |
| 47 | Diff-pair maximum uncoupled length | Altium DiffPairsRouting MAXUNCOUPLEDLENGTH | diff_pair_uncoupled | kind | ODR (500 mil, no pairs) | golden/odrive-v3/src/v3/PCB.PcbDoc record `DiffPairsRouting` |
| 48 | Diff-pair length mismatch | Eagle dpMaxLengthDifference | skew [within_diff_pairs] | kind | EEZ, MOT (10 mm) | golden-rules/*.txt (per-board values; keys resolved in section 6.2) |
| 49 | Impedance target | Altium Width MINIMP, MAXIMP, FAVIMP; impedance-named KiCad netclasses (jetson `Z100_ext` ... `Z50`; names only, no ohm values) | impedance_target | kind | ODR (50 ohm) | golden/odrive-v3/src/v3/PCB.PcbDoc rec Width (MINIMP=50.000000, MAXIMP=50.000000, FAVIMP=50.000000); netclass names only: golden-rules/jetson-nano-baseboard.txt:12-19 |
| 50 | Stackup and materials | KiCad (setup (stackup)); .pro CopperLayerCount, BoardThickness; Eagle layerSetup, mtCopper, mtIsolate | stackup | kind | BMS, EEZ, EVB, GLA, HRF, JET, MOT, SHE, SPD | golden-rules/*.txt (per-board values; keys resolved in section 6.2) |
| 51 | Routing topology | Altium RoutingTopology | routing_topology | kind | ODR (Shortest) | golden-rules/odrive-v3.txt:24 |
| 52 | Routing priority | Altium RoutingPriority | routing_priority | kind | none | golden-rules/odrive-v3.txt:30 |
| 53 | Layers allowed for routing | Altium RoutingLayers | routing_layers | kind | ODR (Mid 1, Mid 2 forbidden) | golden-rules/*.txt (per-board values; keys resolved in section 6.2) |
| 54 | Corner style and angle check | Altium RoutingCorners; Eagle checkAngle | routing_corners | kind | ODR (45°, 100 mil setback) | golden-rules/odrive-v3.txt:28 |
| 55 | Fanout control | Altium FanoutControl | fanout_control | kind | ODR (4 rules) | golden-rules/odrive-v3.txt:16-19 |
| 56 | Drill layer pairs | Altium LayerPairs | via_span.enforce_pairs | kind | ODR | golden-rules/odrive-v3.txt:21 |
| 57 | Testpoint style | Altium FabricationTestpoint, AssemblyTestpoint | testpoint_style | kind | ODR | golden-rules/odrive-v3.txt:4-5 |
| 58 | Testpoint usage | Altium FabricationTestPointUsage, AssemblyTestPointUsage | testpoint_usage | kind | ODR | golden-rules/odrive-v3.txt:6-7 |
| 59 | Unrouted net | Altium UnRoutedNet; KiCad unconnected_items (default error, unchanged) | unrouted_net | kind | ODR | golden-rules/odrive-v3.txt:36 |
| 60 | Short circuit | Altium ShortCircuit ALLOWED=FALSE; KiCad shorting_items (default error, unchanged) | short_circuit | kind | ODR | golden-rules/odrive-v3.txt:37 |
| 61 | Antennae and dangling copper | Altium NetAntennae; KiCad track_dangling, via_dangling (stored warning, the default); Eagle checkWireStubs | dangling_track + dangling_via | kind | ODR (0 mil tolerance), MOT (Eagle 1), KiCad 6+ boards (warning) | golden-rules/odrive-v3.txt:9; golden-rules/bms-c1.txt:4 |
| 62 | Unpoured or modified polygon | Altium UnpouredPolygon | unpoured_zone | kind | ODR | golden-rules/odrive-v3.txt:8 |
| 63 | Keepout (rule) areas | KiCad keepout zones; Eagle tRestrict/bRestrict/vRestrict with checkRestrict | disallow | kind | EEZ, EVB, GLA, HRF, JET, SHE, SPD (KiCad); EEZ (Eagle, checkRestrict 1) | golden-rules/*.txt (per-board values; keys resolved in section 6.2) |
| 64 | Net ties | KiCad net_tie_pad_groups; Altium ComponentClearance scoped to a net-tie footprint | net_tie (Altium: courtyard_clearance relaxed for the net-tie footprint) | kind | BMS 1, GLA (revD) 1, MOT 7; ODR | golden-rules/odrive-v3.txt:23 |
| 65 | Footprint- or component-scoped rule exception | Altium scope queries HasFootprint(...), Not InComponent('J1') | @scope.expr with footprint() / component() and @priority | model | ODR | golden-rules/odrive-v3.txt:13, :23, :38 |
| 66 | DRC severity and check selection | KiCad rule_severities; Eagle check*, maxErrors; Altium RULESETTOCHECK, ONLINERULESETTOCHECK | @severity, @enabled; drc_engine_settings | model | BMS, EVB, GLA, HRF, JET, MOT, SHE, SPD (KiCad); EEZ, MOT (Eagle); ODR | golden-rules/*.txt (per-board values; keys resolved in section 6.2) |
| 67 | ERC severity configuration | KiCad erc.rule_severities; .pro ERC_TestSimilarLabels; Altium [Electrical Rules Check] | @severity of erc_* kinds | model | BMS, EVB, GLA, HRF, JET, MOT, SHE, SPD; EEZ, EVB, GLA, SHE, SPD (.pro flag 1); ODR | golden-rules/*.txt (per-board values; keys resolved in section 6.2) |
| 68 | ERC pin-conflict matrix | KiCad erc.pin_map; Altium [ERC Connection Matrix] | erc_pin_type_conflict.matrix | kind | MOT (edited in 3 copies); ODR | golden-rules/*.txt (per-board values; keys resolved in section 6.2) |
| 69 | Waivers (exclusions) | KiCad erc_exclusions, drc_exclusions | @waiver (finding status dismissed) | model | MOT 27 ERC, SHE 1 ERC | golden-rules/*.txt (per-board values; keys resolved in section 6.2) |
| 70 | Arc approximation tolerance | rules.max_error | drc_engine_settings.arc_tolerance | kind | BMS, EVB, GLA, HRF, JET, MOT, SHE, SPD (0.005 mm) | golden-rules/*.txt (per-board values; keys resolved in section 6.2) |
| 71 | Length measured with stackup height | rules.use_height_for_length_calcs | length.include_via_height | kind | BMS, EVB, GLA, HRF, JET, MOT, SHE, SPD (true) | golden-rules/*.txt (per-board values; keys resolved in section 6.2) |
| 72 | Router presets (widths, vias, pairs) | track_widths, via_dimensions, diff_pair_dimensions; KiCad 5 user_trace_width, user_via; .pro TrackWidthN, ViaDiameterN, dPair*N | router_presets | kind | BMS, EEZ, EVB, GLA, HRF, JET, MOT, SHE | golden-rules/*.txt (per-board values; keys resolved in section 6.2) |
| 73 | Teardrop generator settings | teardrop_options, teardrop_parameters | teardrop_style | kind | BMS, EVB (Rev L), GLA, JET, MOT, SHE | golden-rules/*.txt (per-board values; keys resolved in section 6.2) |
| 74 | Zone fill behaviour | zones_allow_external_fillets | zone_fill_style.external_fillets | kind | none | golden-rules/*.txt (per-board values; keys resolved in section 6.2) |
| 75 | Creepage groove width | rules.min_groove_width | creepage.min_groove_width | kind | none | golden-rules/*.txt (per-board values; keys resolved in section 6.2) |
| 76 | Tuning profiles | netclass tuning_profile; tuning_profiles_impedance_geometric | impedance_target + tuning_profile_consistency | kind | none | golden-rules/glasgow.txt:25, :38 |
| 77 | Eagle pad-shape overrides | psTop, psBottom, psFirst, psElongationLong, psElongationOffset | pad_shape_rule | kind | EEZ, MOT | golden-rules/*.txt (per-board values; keys resolved in section 6.2) |
| 78 | Eagle SMD roundness | srRoundness, srMinRoundness, srMaxRoundness | pad_shape_rule.roundness | kind | none | golden-rules/*.txt (per-board values; keys resolved in section 6.2) |
| 79 | Eagle vias in SMD pads | mnLayersViaInSmd | vias_under_smd (meaning of `mnLayersViaInSmd` unknown) | provisional | EEZ, MOT (2) | golden-rules/*.txt (per-board values; keys resolved in section 6.2) |
| 80 | Eagle optional checks and limits | checkGrid, checkAngle, checkFont, checkRestrict, checkStop, checkValues, checkNames, checkWireStubs, checkPolygonWidth, useDiameter, maxErrors | board_off_grid, routing_corners, text_font, disallow, solder_mask_bridge, part_value_present, erc_annotation, dangling_track, zone_min_width (each @enabled); annular_ratio (`useDiameter`); drc_engine_settings (`maxErrors`) | kind | EEZ, MOT | golden-rules/*.txt (per-board values; keys resolved in section 6.2) |

### 6.2 Native keys in the golden-board extracts

Every distinct rule-bearing key found in the ten `golden-rules/<board>.txt` extracts, resolved through the tables of
section 4 (DRC and ERC severity keys include the ones the extract hides because they are stored as `error`):

| Native family (golden-rules extracts) | Distinct keys | -> canonical kind | -> model field | Excluded (not a rule) | Provisional | Unmapped |
|---|---|---|---|---|---|---|
| Altium RULEKIND | 32 | 32 | 0 | 0 | 0 | 0 |
| Eagle `<param>` names | 73 | 70 | 1 | 0 | 1 | 1 |
| Eagle class fields | 1 | 1 | 0 | 0 | 0 | 0 |
| KiCad 5 (net_class) keys | 8 | 8 | 0 | 0 | 0 | 0 |
| KiCad DRC severity keys | 35 | 35 | 0 | 0 | 0 | 0 |
| KiCad ERC severity keys | 24 | 24 | 0 | 0 | 0 | 0 |
| KiCad board (setup) tokens | 97 | 36 | 0 | 61 | 0 | 0 |
| KiCad netclass fields | 17 | 10 | 2 | 5 | 0 | 0 |
| KiCad project lines | 8 | 7 | 1 | 0 | 0 | 0 |
| KiCad project rules.* | 25 | 25 | 0 | 0 | 0 | 0 |
| **Total** | 320 | 248 | 4 | 66 | 1 | 1 |

Coverage: 252 of 254 rule-bearing native keys map to a canonical kind or model field (99.2 %); 1 provisional, 1 unmapped; 66 further keys are plot, display or creation settings and are excluded by definition (the 61 board `(setup)` tokens are listed below; the other 5 are the netclass display fields of section 4.1.6 rows 13-17). Exceptions:

| Family | Key | Status | Boards |
|---|---|---|---|
| Eagle `<param>` names | `mdSmdStop` | UNMAPPED: none: meaning unknown (Q-E2) | eez-psu-h24005, moteus |
| Eagle `<param>` names | `mnLayersViaInSmd` | provisional | eez-psu-h24005, moteus |

Both exceptions are Eagle parameters whose meaning no source gives (Q-E2); both appear only as constant values in the
EEZ and moteus Eagle files (`mdSmdStop` 0 mil, `mnLayersViaInSmd` 2).

Per board: bms-c1, esp32-evb, glasgow, hackrf-one, jetson-nano-baseboard, sh-esp32 and speeduino-v04 (KiCad only)
map completely. eez-psu-h24005 and moteus map completely on the KiCad side; on the Eagle side everything maps except
`mdSmdStop` (unmapped) and `mnLayersViaInSmd` (provisional). odrive-v3 (Altium) maps completely: 37 of 37 records and
32 of 32 kinds, plus its DRC options and `.PrjPcb` ERC settings, with the per-type ERC levels waiting on Q-A6. The kind-level figure for K1 is therefore
80 of 80 golden constraint kinds mapped (one provisional), and 252 of 254 rule-bearing keys.

KiCad board `(setup)` tokens, all 97 (the union of the `setup keys:` lines of every golden `.kicad_pcb` section; 58
occur in KiCad 4 and 5 files, 67 in newer files, 28 in both). Each is resolved by the row named.

Mapped to a canonical kind (36):

| Token(s) | Resolved by |
|---|---|
| `trace_clearance` | 4.1.15 row 1 |
| `zone_clearance` | 4.1.15 row 2 |
| `trace_min` | 4.1.15 row 4 |
| `user_trace_width` | 4.1.15 row 6 |
| `via_size` | 4.1.15 row 9 |
| `via_drill` | 4.1.15 row 10 |
| `via_min_size` | 4.1.15 row 11 |
| `via_min_drill` | 4.1.15 row 12 |
| `user_via` | 4.1.15 row 13 |
| `uvia_size` | 4.1.15 row 14 |
| `uvia_drill` | 4.1.15 row 15 |
| `uvias_allowed` | 4.1.15 row 16 |
| `uvia_min_size` | 4.1.15 row 17 |
| `uvia_min_drill` | 4.1.15 row 18 |
| `pad_to_mask_clearance` | 4.1.15 row 26 |
| `solder_mask_min_width` | 4.1.15 row 27 |
| `viasonmask` | 4.1.15 row 30 |
| `allow_soldermask_bridges_in_footprints` | 4.1.20 row 1 |
| `tenting` | 4.1.20 row 2 |
| `covering`, `plugging` | 4.1.20 row 3 |
| `capping`, `filling` | 4.1.20 row 4 |
| `front`, `back` | 4.1.20 row 5 |
| `pad_to_paste_clearance_ratio` | 4.1.20 row 6 |
| `stackup` | 4.1.20 row 7 |
| `layer`, `type` | 4.1.20 row 8 |
| `thickness` | 4.1.20 row 9 |
| `material` | 4.1.20 row 10 |
| `epsilon_r` | 4.1.20 row 11 |
| `loss_tangent` | 4.1.20 row 12 |
| `color` | 4.1.20 row 13 |
| `copper_finish` | 4.1.20 row 14 |
| `dielectric_constraints` | 4.1.20 row 15 |

Excluded, not a rule (61: plot, display, creation and coordinate settings):

| Token(s) | Resolved by |
|---|---|
| `zone_45_only` | 4.1.15 row 3 |
| `last_trace_width` | 4.1.15 row 5 |
| `segment_width` | 4.1.15 row 7 |
| `edge_width` | 4.1.15 row 8 |
| `pcb_text_width` | 4.1.15 row 19 |
| `pcb_text_size` | 4.1.15 row 20 |
| `mod_edge_width` | 4.1.15 row 21 |
| `mod_text_width` | 4.1.15 row 22 |
| `mod_text_size` | 4.1.15 row 23 |
| `pad_size` | 4.1.15 row 24 |
| `pad_drill` | 4.1.15 row 25 |
| `aux_axis_origin`, `grid_origin` | 4.1.15 row 28 |
| `pcbplotparams`, `creategerberjobfile`, `excludeedgelayer`, `hpglpendiameter`, `hpglpennumber`, `hpglpenoverlay`, `hpglpenspeed`, `layerselection`, `linewidth`, `mirror`, `mode`, `outputdirectory`, `outputformat`, `padsonsilk`, `plotframeref`, `plotinvisibletext`, `plotreference`, `plotvalue`, `psa4output`, `psnegative`, `scaleselection`, `useauxorigin`, `usegerberadvancedattributes`, `usegerberattributes`, `usegerberextensions` | 4.1.15 row 29 |
| `subtractmaskfromsilk` | 4.1.15 row 31 |
| `drillshape` | 4.1.15 row 32 |
| `visible_elements` | 4.1.15 row 33 |
| `crossoutdnponfab`, `dashed_line_dash_ratio`, `dashed_line_gap_ratio`, `disableapertmacros`, `dxfimperialunits`, `dxfpolygonmode`, `dxfusepcbnewfont`, `hidednponfab`, `pdf_back_fp_property_popups`, `pdf_front_fp_property_popups`, `pdf_metadata`, `pdf_single_document`, `plot_black_and_white`, `plot_on_all_layers_selection`, `plotfptext`, `plotpadnumbers`, `sketchdnponfab`, `sketchpadsonfab`, `svgprecision`, `svguseinch` | 4.1.20 row 16 |

### 6.3 odrive-v3 Altium rules (parsed from `Rules6`, all 37 records, 32 kinds)

| # | NAME | RULEKIND | Scope 1 / scope 2 | PRIORITY | Canonical kind(s) | Canonical scope | Provenance |
|---|---|---|---|---|---|---|---|
| 1 | MinimumAnnularRing | MinimumAnnularRing | `All` / `All` | 1 | annular_width | expr: `all` | golden-rules/odrive-v3.txt:3 |
| 2 | FabricationTestpoint | FabricationTestpoint | `All` / `All` | 1 | testpoint_style | expr: `all` | golden-rules/odrive-v3.txt:4 |
| 3 | AssemblyTestpoint | AssemblyTestpoint | `All` / `All` | 1 | testpoint_style | expr: `all` | golden-rules/odrive-v3.txt:5 |
| 4 | AssemblyTestPointUsage | AssemblyTestPointUsage | `All` / `All` | 1 | testpoint_usage | expr: `all` | golden-rules/odrive-v3.txt:6 |
| 5 | FabricationTestPointUsage | FabricationTestPointUsage | `All` / `All` | 1 | testpoint_usage | expr: `all` | golden-rules/odrive-v3.txt:7 |
| 6 | UnpouredPolygon | UnpouredPolygon | `All` / `All` | 1 | unpoured_zone | expr: `all` | golden-rules/odrive-v3.txt:8 |
| 7 | NetAntennae | NetAntennae | `All` / `All` | 1 | dangling_track | expr: `all` | golden-rules/odrive-v3.txt:9 |
| 8 | SilkToBoardRegionClearance | SilkToBoardRegionClearance | `All` / `All` | 1 | silk_to_edge_clearance | expr: `all` | golden-rules/odrive-v3.txt:10 |
| 9 | SilkToSilkClearance | SilkToSilkClearance | `All` / `All` | 1 | silk_clearance | expr: `all` | golden-rules/odrive-v3.txt:11 |
| 10 | SilkToSolderMaskClearance | SilkToSolderMaskClearance | `IsPad` / `All` | 1 | silk_to_mask_clearance | expr: `type(pad)` | golden-rules/odrive-v3.txt:12 |
| 11 | MinimumSolderMaskSliver | MinimumSolderMaskSliver | `InPadClass('All Pads') And Not InComponent('J1')` / `All` | 1 | solder_mask_sliver | expr: `pad_class('All Pads') and not component('J1')` | golden-rules/odrive-v3.txt:13 |
| 12 | HoleToHoleClearance | HoleToHoleClearance | `All` / `All` | 1 | hole_to_hole + microvia_stack | expr: `all` | golden-rules/odrive-v3.txt:14 |
| 13 | DiffPairsRouting | DiffPairsRouting | `All` / `All` | 1 | diff_pair_gap + diff_pair_uncoupled + track_width + impedance_target | expr: `all` | golden-rules/odrive-v3.txt:15 |
| 14 | Fanout_LCC | FanoutControl | `IsLCC` / `All` | 1 | fanout_control | expr: `package_style('LCC')` | golden-rules/odrive-v3.txt:16 |
| 15 | Fanout_SOIC | FanoutControl | `IsSOIC` / `All` | 2 | fanout_control | expr: `package_style('SOIC')` | golden-rules/odrive-v3.txt:17 |
| 16 | Fanout_Small | FanoutControl | `(CompPinCount < 5)` / `All` | 3 | fanout_control | expr: `(pin_count() < 5)` | golden-rules/odrive-v3.txt:18 |
| 17 | Fanout_Default | FanoutControl | `All` / `All` | 4 | fanout_control | expr: `all` | golden-rules/odrive-v3.txt:19 |
| 18 | HoleSize | HoleSize | `All` / `All` | 1 | hole_size | expr: `all` | golden-rules/odrive-v3.txt:20 |
| 19 | LayerPairs | LayerPairs | `All` / `All` | 1 | via_span | expr: `all` | golden-rules/odrive-v3.txt:21 |
| 20 | Height | Height | `All` / `All` | 1 | component_height | expr: `all` | golden-rules/odrive-v3.txt:22 |
| 21 | ComponentClearance | ComponentClearance | `HasFootprint('NetTie_0.25mm')` / `All` | 1 | courtyard_clearance | expr: `footprint('NetTie_0.25mm')` | golden-rules/odrive-v3.txt:23 |
| 22 | RoutingTopology | RoutingTopology | `All` / `All` | 1 | routing_topology | expr: `all` | golden-rules/odrive-v3.txt:24 |
| 23 | SolderMaskExpansion | SolderMaskExpansion | `All` / `All` | 1 | solder_mask_expansion + via_protection | expr: `all` | golden-rules/odrive-v3.txt:25 |
| 24 | PasteMaskExpansion | PasteMaskExpansion | `All` / `All` | 1 | solder_paste_margin | expr: `all` | golden-rules/odrive-v3.txt:26 |
| 25 | RoutingVias | RoutingVias | `All` / `All` | 1 | via_diameter + hole_size + via_span | expr: `all` | golden-rules/odrive-v3.txt:27 |
| 26 | RoutingCorners | RoutingCorners | `All` / `All` | 1 | routing_corners | expr: `all` | golden-rules/odrive-v3.txt:28 |
| 27 | RoutingLayers | RoutingLayers | `All` / `All` | 1 | routing_layers | expr: `all` | golden-rules/odrive-v3.txt:29 |
| 28 | RoutingPriority | RoutingPriority | `All` / `All` | 1 | routing_priority | expr: `all` | golden-rules/odrive-v3.txt:30 |
| 29 | Width | Width | `All` / `All` | 1 | track_width + impedance_target | expr: `all` | golden-rules/odrive-v3.txt:31 |
| 30 | PolygonConnect | PolygonConnect | `All` / `All` | 1 | zone_connection + thermal_relief_gap + thermal_spoke_width + thermal_spoke_count | expr: `all` | golden-rules/odrive-v3.txt:32 |
| 31 | PlaneClearance | PlaneClearance | `All` / `All` | 1 | plane_clearance | expr: `all` | golden-rules/odrive-v3.txt:33 |
| 32 | PlaneConnect | PlaneConnect | `All` / `All` | 1 | zone_connection + thermal_relief_gap + thermal_spoke_width + thermal_spoke_count | expr: `all` | golden-rules/odrive-v3.txt:34 |
| 33 | Clearance | Clearance | `All` / `All` | 1 | clearance + hole_clearance + physical_clearance | expr: `all`; pair: [different_net] | golden-rules/odrive-v3.txt:35 |
| 34 | UnRoutedNet | UnRoutedNet | `All` / `All` | 1 | unrouted_net | expr: `all` | golden-rules/odrive-v3.txt:36 |
| 35 | ShortCircuit | ShortCircuit | `All` / `All` | 1 | short_circuit + tracks_crossing | expr: `all` | golden-rules/odrive-v3.txt:37 |
| 36 | ComponentClearance_1 | ComponentClearance | `HasFootprint('1704936') or HasFootprint('1868076')` / `All` | 2 | courtyard_clearance | expr: `footprint('1704936') or footprint('1868076')` | golden-rules/odrive-v3.txt:38 |
| 37 | ComponentClearance_2 | ComponentClearance | `All` / `All` | 3 | courtyard_clearance | expr: `all` | golden-rules/odrive-v3.txt:39 |

Besides the 37 records, ODrive's DRC options (section 4.2.4) and `.PrjPcb` ERC levels, connection matrix and class
generation flags (section 4.2.5) map as listed there; ODrive has no differential pairs (`DifferentialPairs6/Data` holds
0 records), so its DiffPairsRouting rule constrains nothing.

### 6.4 What the golden boards say about the taxonomy

- No golden board has a `.kicad_dru` file (search of `golden/*/src` found none), so the DRU language, explicit
  priority and binary scopes are exercised only by ODrive's Altium rules and the Eagle class matrices. K1 needs
  synthetic rule-probe boards (docs/13-eda-knowledge-ingestion.md) for the rest.
- High-voltage intent exists only as netclass clearance (EEZ `Mains` 1.5 mm, golden-rules/eez-psu-h24005.txt:38;
  bms-c1 `HV` 0.3 mm, golden-rules/bms-c1.txt:12); no creepage rule, `min_groove_width` 0 mm. `creepage` is in the
  taxonomy but untested on real data.
- Impedance intent exists only in netclass names and widths (jetson `Z100`, `Z90`, `Z85`, `Z50`,
  golden-rules/jetson-nano-baseboard.txt:12-19) and in ODrive's Width impedance fields (50 ohm,
  golden/odrive-v3/src/v3/PCB.PcbDoc rec Width): `impedance_target` is needed, and a class-name heuristic is not a
  rule (it would be a T3 suggestion).
- Two boards set minimum via diameters that KiCad's validator calls ineffective (bms-c1 0.2 mm against 0.4 mm;
  moteus 0.381 mm against 0.4572 mm; golden-rules/bms-c1.txt:3; golden-rules/moteus.txt:3): the kernel should compute
  the effective floor `max(via_diameter.min, hole_size.min + 2 x annular_width.min)` and report the difference.
- esp32-evb gives up to 169 pads a local clearance of 0.0508 mm against a 0.1778 mm netclass clearance. By the code
  path in P9 the pad override wins over the netclass and is floored only at `min_clearance` (0 mm on that board), so the
  effective clearance on those pads is 0.0508 mm (inferred from pcbnew/drc/drc_engine.cpp:1180-1247 and
  pcbnew/pad.cpp:1913-1921†; confirm in L1). This settles the question left open in the golden-usage notes.
- Checks are often switched off (hackrf-one edge clearance 0 mm with `invalid_outline` ignored; sh-esp32 all silk and
  text checks ignored; jetson ERC `pin_to_pin` and `power_pin_not_driven` ignored; seven boards ignore holes inside
  courtyards; golden-rules/hackrf-one.txt:3-4; golden-rules/sh-esp32.txt:17; golden-rules/jetson-nano-baseboard.txt:21):
  `severity: ignore` and `enabled` must survive import, and the kernel should report disabled checks.

## 7. Plug-in points

### 7.1 Fab capability floors (W7)

Each `fabs/<fab>.yaml` entry names a canonical kind and the bound the fab can make, per process. The loader turns
each entry into a rule with `origin: fab_floor`, `scope.fab` and `scope.process` set, `provenance` = the capability
page URL and `verified_at`, and applies it by P8 (a floor never relaxes a design rule; only a recorded waiver may go
below it). Kinds a fab file is expected to fill: `track_width`, `clearance`, `hole_size`, `via_diameter`,
`annular_width`, `hole_to_hole`, `hole_clearance`, `edge_clearance`, `hole_aspect_ratio`, `via_span`,
`microvia_stack`, `solder_mask_sliver`, `solder_mask_expansion`, `silk_line_width`, `text_height`, `text_thickness`,
`board_size`, `stackup`, `vias_under_smd` (via-in-pad), `copper_sliver`. Shape only, no values:

```yaml
# fabs/<fab>.yaml (W7) - shape of one entry; values are placeholders
fab: <fab_id>
process: <process_id>
verified_at: <YYYY-MM-DD>
source_url: <capability page>
floors:
  - kind: track_width
    layers: [outer]
    constraint: {min: <value>, unit: mm}
  - kind: hole_size
    constraint: {min: <value>, max: <value>, unit: mm}
```

### 7.2 Interface rules (W8)

W8 rules are ordinary `rule.schema.json` records. They attach to canonical kinds through `kind` and to objects
through the existing selectors (`interface`, `part`, `package`, `stage`) plus `scope.expr` for nets or components the
netlist classifier tagged, with `origin: explicit`. A W8 rule never relaxes a floor (P8) and never overrides an
item override (P3). Shape only, no rule content; angle-bracket placeholders stand for numbers or strings:

```json
{
  "id": "if.<interface>.<rule>",
  "kind": "skew",
  "scope": {"interface": "usb2", "stage": "routing", "expr": "diffpair('<base net>')"},
  "constraint": {"max": "<value>", "unit": "ps", "options": {"within_diff_pairs": true}},
  "origin": "explicit",
  "priority": "<n>",
  "severity": "error",
  "title": "<title>",
  "check": "pcbbrain_checks.drc.skew",
  "provenance": [{"doc": "<app note>", "section": "<section>", "url": "<url>"}],
  "confidence": "<high|medium|low>"
}
```

Rules that are not geometric (for example "a 22 ohm series resistor on each USB data line") keep binding to a check
callable only, without `kind`, as `rule.schema.json` already allows.

## 8. Open questions and things to verify locally (L1)

### 8.1 Altium (verify against the local Altium install, by scripting, never by inspecting binaries)

- **Q-A1** altium.com was blocked: every AD-keyed fact (category lists, UI names, units) comes from search summaries.
  In L1, list every rule kind with its UI name, category and parameter names from the PCB Rules and Constraints Editor
  through Altium's scripting interface (an `IterateRules`-style DelphiScript, W1), and compare with section 4.2.1.
  UI names currently unverified (GK): Solder Mask Expansion, Paste Mask Expansion, Power Plane Connect Style, Power
  Plane Clearance, Polygon Connect Style, Fabrication / Assembly Testpoint Style and Usage, Silk To Board Region
  Clearance, Routing Neck-Down, Wire Bonding, Layer Stack.
- **Q-A2** Parameter keys were found in no record for 26 kinds: UnConnectedPin, Creepage, RoutingNeckDown, SMDToPlane,
  SMDNeckDown, AcuteAngle, ParallelSegment, Length, StubLength, MaximumViaCount, BackDrilling, ReturnPath,
  ComponentOrientations, PermittedLayers, NetsToIgnore, the 12 signal-integrity kinds other than SupplyNets,
  LayerStack, WireBonding, ZAxisClearance. Create one rule of each in a probe PcbDoc and read the `Rules6` record.
- **Q-A3** Unary or binary is sourced only for Clearance and Return Path; check the two-scope editor for every kind.
- **Q-A4** Meanings of `COLLISIONCHECKMODE=3` (ComponentClearance), `SIDE=3` (FabricationTestpoint), `VALID=0`
  (testpoint usage); units of `DEADCOPPERMINAREA` and of the suffix-less SMDToCorner `DISTANCE`.
- **Q-A5** AD1 says 65 rule types; there are 69 real `TRuleKind` values. Reconcile (likely LayerStack, WireBonding,
  ZAxisClearance and one more are outside the count).
- **Q-A6** Map `[Electrical Rules Check]` `Type1`-`Type111` to violation names, and confirm the row and column labels of
  the 17 x 17 `[ERC Connection Matrix]`; until then every Altium ERC mapping is lossy.
- **Q-A7** Does Altium still apply the legacy scope-kind precedence (`TScopeKind`) next to `PRIORITY`?
- **Q-A8** Categories of SilkToBoardRegionClearance, RoutingNeckDown, WireBonding, ZAxisClearance and LayerStack.
- **Q-A9** Query semantics: is a binary rule matched in both orders; do `And` and `Or` follow the usual precedence;
  what do `Like`, `Div`, `Mod`, `Between`, bare net names and `OnMultiLayer` mean exactly.
- **Q-A10** Keepouts: confirm they are primitives with a `KEEPOUT` flag and not a rule kind; confirm how Altium
  handles net-tie components (the canonical `net_tie` has no Altium mapping yet).
- **Q-A11** Decode the ODrive layer stack (`Board6/Data`) so per-layer Width, DiffPairsRouting and RoutingLayers keys
  map to `cu:<k>`.
- **Q-A12** Project rule D6 forbids decompiling or instrumenting Altium binaries. The draft cited three third-party
  sources derived from decompiled Altium .NET assemblies (an AD26 type dump and a crate built on it). Those citations
  were removed in the same session (docs/11 Q26, default: do not cite such sources). Facts that rested only on them
  are marked "unverified (L1: confirm via Altium scripting, Q-A12)": confirm or correct each one in L1 with Altium's
  own scripting interface (TRuleKind enumeration and rule-object properties via DelphiScript, docs/13 § 13.2b),
  which is the D6-approved route. File-format observations (KiCad importer, golden Rules6 records, public ASCII
  PcbDoc and Rules6 dumps, Altium's own DelphiScript sample W1) remain cited.

### 8.2 KiCad (verify with `kicad-cli` and probe boards)

- **Q-K1** How to write a non-zero `solder_paste_rel_margin` in a `.kicad_dru`: the parser treats it as a length, the
  pad code divides the stored value by 1000 (pcbnew/drc/drc_rule_parser.cpp:612-619; pcbnew/pad.cpp:2064†).
- **Q-K2** `solder_mask_sliver` is parsed and offered by the graphical editor but no consumer was found in
  `pcbnew/` or `common/`; plotting code elsewhere was not searched.
- **Q-K3** First match (help) or per-value last match (code)? Probe with two rules setting different fields.
- **Q-K4** `rules.min_via_annular_width` default: 0.1 mm (constructor) or 0.5 mm (`PARAM` default),
  pcbnew/board_design_settings.cpp:155, :304-306.
- **Q-K5** The implicit courtyard rule sets its min on the hole-to-hole constraint object
  (pcbnew/drc/drc_engine.cpp:233-236), so the effective default courtyard clearance depends on `MINOPTMAX` defaults
  (include/core/minoptmax.h, not in the checkout).
- **Q-K6** Operator precedence (`||` tighter than `&&`, inferred from grammar.lemon); test `A || B && C`. Also confirm
  that `!` binds tighter than `==` (grammar.lemon:31-32): `!A.NetName == 'GND'` should parse as `(!A.NetName) == 'GND'`.
- **Q-K7** `netclass_patterns` matcher syntax (wildcard or regex), defined outside the checkout.
- **Q-K8** Consumer of `rules.use_height_for_length_calcs`.
- **Q-K9** Whether project migration (project_file.cpp, not in the checkout) consumes `allow_microvias`,
  `allow_blind_buried_vias` and the ERC key `global_label_dangling`; without it they are ignored.
- **Q-K10** Explicit rules escape the board-minimum floor (inferred). Import moteus `hw/controller/controller.brd` and
  check that CAN-to-CAN resolves to 0.127 mm and default pairs to 0.2 mm.
- **Q-K11** Smaller items: `(time_domain yes)` appears in help but not in the parser; `(min -1)` without units would be
  "Missing units"; `DIFF_PAIR_INTRA_SKEW_CONSTRAINT` has no consumer; `DRCE_GENERIC_WARNING` has error severity in the
  map; `DRCE_ZONES_INTERSECT` is dead.
- **Q-K12** Meanings of KiCad 5 `(setup)`, `(net_class)` and `.pro` keys (legacy parsers not in the checkout); the
  mapping of file-format versions to releases; defaults of the KiCad versions that saved each golden board.
- **Q-K13** Constants outside the checkout: `ARC_HIGH_DEF`, zone defaults, advanced-config defaults (DRC epsilon,
  hole wall thickness, sliver thresholds, extra clearance, minimum parallel angle), `MAXIMUM_CLEARANCE`.
- **Q-K14** Release vocabularies exist only for 9.0 and 10.0; versions 6.x-8.x cannot be gated from these sources.

### 8.3 Eagle (verify against a local EAGLE or Fusion Electronics install)

- **Q-E1** No official Eagle documentation was read; every meaning rests on OMX or is unverified (GK-E).
- **Q-E2** Unknown: `mnLayersViaInSmd` (provisionally `vias_under_smd`), `mdSmdStop` (unmapped), the exact scope of
  `mdViaViaSameLayer`, the `useDiameter` bit layout, the base of `dpGapFactor`.
- **Q-E3** `mdDrill`: OMX gives "minimum distance between drill holes" (mapped to `hole_to_hole`), while the
  golden-usage summary groups it with copper-to-hole clearance. Read the Distance tab of the Eagle DRC dialog.
- **Q-E4** Whether Eagle writes `psTop`/`psBottom`/`psFirst` with the codes KiCad assumes (-1 library, 0 square,
  1 round, 2 octagon, 3 long, 4 offset).
- **Q-E5** How Eagle combines the class matrix with the `md*` values (larger wins?), and what a missing matrix entry or
  a 0 mil `md*` means. Canonical default: `md*` are floors, matrix entries are floored explicit rules.
- **Q-E6** The ` *` suffix on rule-set names and the `mtIsolate` indexing are data-driven hypotheses.
- **Q-E7** Eagle schematic net classes: the eeschema Eagle importer is not in the checkout.
- **Q-E8** Fix the golden extractor so unnamed `<designrules>` blocks are captured (session-scratchpad script
  `scripts/golden_rules.py`, line 54, which requires a `name=` attribute).

### 8.4 General

- **Q-G1** DesignIR (`schemas/design-graph.md`) lacks pads, holes, text, graphics, courtyards, rule areas, net
  chains, component classes, groups and layer ids as typed objects; W6 must add the types marked (+) in section 2.3
  before scope expressions can be evaluated.
- **Q-G2** No golden board exercises the DRU language; build synthetic probe boards (one kind per board, parameters
  swept) for the KiCad oracle runs.
- **Q-G3** `silk_line_width` and `board_size` are project-defined for W7; confirm against the fabs' capability pages.
- **Q-G4** Altium PCB rules carry no severity; confirm how Altium reports each kind (all violations, or per-kind
  levels) so `inherit` resolves correctly.
- **Q-G5** `floor_exempt` on imported KiCad DRU rules reproduces KiCad but lets explicit rules sit below the board
  minimum; the kernel reports them. If the owner prefers floors to always apply, drop the flag on import (one-line
  change in the importer; record in docs/11).

## Appendix A. KiCad expression properties

All 343 property registrations at this commit, by registering class (derived classes inherit unless they
mask a property; masks are listed in pcbnew/*.cpp† property-registration blocks). Every name maps to
`prop('<Display Name>')`; the last column lists the names with a first-class canonical form. Registered but not
readable from an expression: `PAD.Top_Post-machining`, `PAD.Top_Post-machining_Size`, `PAD.Bottom_Post-machining`, `PAD.Bottom_Post-machining_Size`, `PAD.Bottom_Backdrill_Must-Cut`, `PAD.Top_Backdrill_Must-Cut`, `PCB_DRILL_CHART.Non-plated_Holes`, `PCB_TUNING_PATTERN.Corner_Radius_%`, `PCB_TUNING_PATTERN.Single-sided`, `PCB_VIA.Bottom_Backdrill_Must-Cut`, `PCB_VIA.Top_Backdrill_Must-Cut`, `PCB_VIA.Front_Post-machining`, `PCB_VIA.Front_Post-machining_Size`, `PCB_VIA.Front_Post-machining_Depth`, `PCB_VIA.Front_Post-machining_Angle`, `PCB_VIA.Back_Post-machining`, `PCB_VIA.Back_Post-machining_Size`, `PCB_VIA.Back_Post-machining_Depth`, `PCB_VIA.Back_Post-machining_Angle`, `PCB_VIA_STACK.Copper-filled`, `PCB_VIA_STITCH.Seed`.

| Registering class | Count | Expression names (A.&lt;name>) | First-class canonical form (others: prop('&lt;Display Name>')) | Provenance |
|---|---|---|---|---|
| `BOARD_CONNECTED_ITEM` | 14 | Layer, Net, Net_Class, NetClass, NetName, Enable_Teardrops, Best_Length_Ratio, Max_Length, Best_Width_Ratio, Max_Width, Curved_Teardrops, Prefer_Zone_Connections, Allow_Teardrops_To_Span_Two_Tracks, Max_Width_Ratio | Layer -> on_layer(); Net -> net(); Net_Class -> netclass(); NetClass -> netclass(); NetName -> net() | pcbnew/board_connected_item.cpp:223† (first registration) |
| `BOARD_ITEM` | 5 | Parent, Position_X, Position_Y, Layer, Locked | Layer -> on_layer() | pcbnew/board_item.cpp:540† (first registration) |
| `EDA_GROUP` | 1 | Name | Name -> group() / prop('Name') | pcbnew/pcb_group.cpp:631† (first registration) |
| `EDA_SHAPE` | 30 | Shape, Start_X, Start_Y, Center_X, Center_Y, Radius, End_X, End_Y, Width, Height, Corner_Radius, Major_Radius, Minor_Radius, Ellipse_Rotation, Arc_Start_Angle, Arc_End_Angle, Line_Width, Line_Style, Line_Color, Angle, Fill, Fill_Color, Start_Shape, Start_Length, Start_Width, Start_Stroke_Width, End_Shape, End_Length, End_Width, End_Stroke_Width | Width -> track_width (value) on tracks; Height -> component_height (value) | common/eda_shape.cpp:4746† (first registration) |
| `EDA_TEXT` | 15 | Orientation, Text, Font, Auto_Thickness, Thickness, Italic, Bold, Mirrored, Visible, Width, Height, Horizontal_Justification, Vertical_Justification, Color, Hyperlink | Orientation -> component_orientation (value); Width -> track_width (value) on tracks; Height -> component_height (value) | common/eda_text.cpp:1376† (first registration) |
| `FOOTPRINT` | 21 | Layer, Orientation, Scale_X, Scale_Y, Reference, Library_Link, Library_Description, Keywords, Height, Component_Class, Footprint_Type, Not_in_Schematic, Exclude_From_Position_Files, Exclude_From_Bill_of_Materials, Exclude_From_Simulation, Do_not_Populate, Exempt_From_Courtyard_Requirement, Clearance_Override, Solderpaste_Margin_Override, Solderpaste_Margin_Ratio_Override, Zone_Connection_Style | Layer -> on_layer(); Orientation -> component_orientation (value); Reference -> component(); Library_Link -> footprint(); Height -> component_height (value); Component_Class -> component_class(); Exempt_From_Courtyard_Requirement -> courtyard_required exemption; Clearance_Override -> clearance [item_override]; Solderpaste_Margin_Override -> solder_paste_margin.abs [item_override]; Solderpaste_Margin_Ratio_Override -> solder_paste_margin.ratio [item_override]; Zone_Connection_Style -> zone_connection [item_default] (checked after rules, P9) | pcbnew/footprint.cpp:5275† (first registration) |
| `PAD` | 39 | Orientation, Pad_Type, Pad_Shape, Pad_Number, Pin_Name, Pin_Type, Simulation_Electrical_Type, Size_X, Size_Y, Corner_Radius_Ratio, Corner_Radius_Size, Hole_Shape, Hole_Size_X, Hole_Size_Y, Top_Post-machining, Top_Post-machining_Size, Top_Counterbore_Depth, Top_Countersink_Angle, Bottom_Post-machining, Bottom_Post-machining_Size, Bottom_Counterbore_Depth, Bottom_Countersink_Angle, Backdrill_Mode, Bottom_Backdrill_Size, Bottom_Backdrill_Must-Cut, Top_Backdrill_Size, Top_Backdrill_Must-Cut, Fabrication_Property, Copper_Layers, Pad_To_Die_Length, Pad_To_Die_Delay, Clearance_Override, Soldermask_Margin_Override, Solderpaste_Margin_Override, Solderpaste_Margin_Ratio_Override, Zone_Connection_Style, Thermal_Relief_Spoke_Width, Thermal_Relief_Spoke_Angle, Thermal_Relief_Gap | Orientation -> component_orientation (value); Pad_Type -> type(pad.&lt;t>); Hole_Size_X -> hole_size (value); Clearance_Override -> clearance [item_override]; Soldermask_Margin_Override -> solder_mask_expansion [item_override]; Solderpaste_Margin_Override -> solder_paste_margin.abs [item_override]; Solderpaste_Margin_Ratio_Override -> solder_paste_margin.ratio [item_override]; Zone_Connection_Style -> zone_connection [item_override]; Thermal_Relief_Spoke_Width -> thermal_spoke_width [item_override]; Thermal_Relief_Spoke_Angle -> thermal_spoke_count.angle [item]; Thermal_Relief_Gap -> thermal_relief_gap [item_override] | pcbnew/pad.cpp:3610† (first registration) |
| `PCB_BARCODE` | 11 | Text, Show_Text, Text_Size, Width, Height, Orientation, Barcode_Type, Error_Correction, Knockout, Margin_X, Margin_Y | Width -> track_width (value) on tracks; Height -> component_height (value); Orientation -> component_orientation (value) | pcbnew/pcb_barcode.cpp:1048† (first registration) |
| `PCB_CONSTRAINT` | 1 | Driving | all via prop() | pcbnew/constraints/pcb_constraint.cpp:497† (first registration) |
| `PCB_DIMENSION_BASE` | 17 | Prefix, Suffix, Override_Text, Value_Mode, Value, Text, Units, Units_Format, Precision, Suppress_Trailing_Zeroes, Arrow_Direction, Line_Thickness, Arrow_Length, Extension_Line_Offset, Text_Position, Keep_Aligned_with_Dimension, Orientation | Orientation -> component_orientation (value) | pcbnew/pcb_dimension.cpp:2154† (first registration) |
| `PCB_DIM_ALIGNED` | 2 | Crossbar_Height, Extension_Line_Overshoot | all via prop() | pcbnew/pcb_dimension.cpp:2351† (first registration) |
| `PCB_DIM_LEADER` | 1 | Text_Frame | all via prop() | pcbnew/pcb_dimension.cpp:2469† (first registration) |
| `PCB_DIM_RADIAL` | 1 | Leader_Length | all via prop() | pcbnew/pcb_dimension.cpp:2426† (first registration) |
| `PCB_DRILL_CHART` | 7 | Show_Totals, Plated_Holes, Non-plated_Holes, Vias, Slots, Backdrills, Castellated | all via prop() | pcbnew/pcb_drill_chart.cpp:512† (first registration) |
| `PCB_DRILL_MAP` | 6 | Offset_X, Offset_Y, Symbol_Size, Hole_Span, Outline_Slots, Guide_Cross | all via prop() | pcbnew/pcb_drill_map.cpp:397† (first registration) |
| `PCB_FIELD` | 1 | Name | Name -> group() / prop('Name') | pcbnew/pcb_field.cpp:383† (first registration) |
| `PCB_GENERATED_TABLE` | 7 | Units, Decimal_Places, Row_Separators, Column_Separators, Row_Height, Font_Name, Font_Size | all via prop() | pcbnew/pcb_generated_table.cpp:1117† (first registration) |
| `PCB_GENERATOR` | 1 | Update_Order | all via prop() | pcbnew/pcb_generator.cpp:218† (first registration) |
| `PCB_GRID_ITEM` | 15 | Grid_Type, Orientation, Extent_X, Extent_Y, Spacing_X, Spacing_Y, Radius_Extent, Radius_Spacing, Phi_Extent, Phi_Spacing, Priority, Tick_Interval, Affects_Cursor_Snap, Affects_Routing, Affects_Placement | Orientation -> component_orientation (value) | pcbnew/pcb_grid_item.cpp:542† (first registration) |
| `PCB_POINT` | 1 | Size | all via prop() | pcbnew/pcb_point.cpp:362† (first registration) |
| `PCB_REFERENCE_IMAGE` | 8 | Associated_Layer, Scale, Transform_Offset_X, Transform_Offset_Y, Orientation, Opacity, Width, Height | Orientation -> component_orientation (value); Width -> track_width (value) on tracks; Height -> component_height (value) | pcbnew/pcb_reference_image.cpp:464† (first registration) |
| `PCB_SHAPE` | 5 | Layer, Number_Box, Thermal_Spoke_Template, Soldermask, Soldermask_Margin_Override | Layer -> on_layer(); Soldermask_Margin_Override -> solder_mask_expansion [item_override] | pcbnew/pcb_shape.cpp:2057† (first registration) |
| `PCB_TABLE` | 12 | Start_X, Start_Y, External_Border, Header_Border, Border_Width, Border_Style, Border_Color, Row_Separators, Cell_Separators, Separators_Width, Separators_Style, Separators_Color | all via prop() | pcbnew/pcb_table.cpp:1109† (first registration) |
| `PCB_TABLECELL` | 4 | Font_Name, Font_Size, Column_Width, Row_Height | all via prop() | pcbnew/pcb_tablecell.cpp:433† (first registration) |
| `PCB_TARGET` | 3 | Size, Width, Shape | Width -> track_width (value) on tracks | pcbnew/pcb_target.cpp:245† (first registration) |
| `PCB_TEXT` | 2 | Knockout, Keep_Upright | all via prop() | pcbnew/pcb_text.cpp:928† (first registration) |
| `PCB_TEXTBOX` | 8 | Knockout, Border, Border_Style, Border_Width, Margin_Left, Margin_Top, Margin_Right, Margin_Bottom | all via prop() | pcbnew/pcb_textbox.cpp:1107† (first registration) |
| `PCB_TRACK` | 7 | Width, Start_X, Start_Y, End_X, End_Y, Soldermask, Soldermask_Margin_Override | Width -> track_width (value) on tracks; Soldermask_Margin_Override -> solder_mask_expansion [item_override] | pcbnew/pcb_track.cpp:3613† (first registration) |
| `PCB_TUNING_PATTERN` | 20 | Layer, Width, Net, Soldermask, Soldermask_Margin_Override, End_X, End_Y, Tuning_Mode, Min_Amplitude, Max_Amplitude, Initial_Side, Min_Spacing, Corner_Radius_%, Target_Length, Target_Delay, Target_Skew, Target_Skew_Delay, Override_Custom_Rules, Single-sided, Rounded | Layer -> on_layer(); Width -> track_width (value) on tracks; Net -> net(); Soldermask_Margin_Override -> solder_mask_expansion [item_override] | pcbnew/generators/pcb_tuning_pattern.cpp:3336† (first registration) |
| `PCB_VIA` | 27 | Automatically_Update_Net, Diameter, Hole, Layer_Top, Layer_Bottom, Via_Type, Front_Tenting, Back_Tenting, Front_Covering, Back_Covering, Front_Plugging, Back_Plugging, Capping, Filling, Backdrill_Mode, Bottom_Backdrill_Size, Bottom_Backdrill_Must-Cut, Top_Backdrill_Size, Top_Backdrill_Must-Cut, Front_Post-machining, Front_Post-machining_Size, Front_Post-machining_Depth, Front_Post-machining_Angle, Back_Post-machining, Back_Post-machining_Size, Back_Post-machining_Depth, Back_Post-machining_Angle | Diameter -> via_diameter (value); Hole -> hole_size (value); Via_Type -> type(via.&lt;t>); Front_Tenting -> via_protection.tent_top [item]; Back_Tenting -> via_protection.tent_bottom [item]; Filling -> via_protection.fill [item] | pcbnew/pcb_track.cpp:3664† (first registration) |
| `PCB_VIA_STACK` | 10 | Start_Layer, End_Layer, Style, Net, Use_Netclass_Values, Via_Diameter, Via_Hole, Pitch, Copper-filled, Capped | Net -> net() | pcbnew/generators/pcb_via_stack.cpp:1645† (first registration) |
| `PCB_VIA_STITCH` | 8 | Size, Drill, Pitch, Layout, Mode, Seed, Net, Guarded_Net | Net -> net() | pcbnew/generators/pcb_via_stitch.cpp:1710† (first registration) |
| `ZONE` | 33 | Position_X, Position_Y, Layer, Priority, Name, Keep_Out_Tracks, Keep_Out_Vias, Keep_Out_Pads, Keep_Out_Zone_Fills, Keep_Out_Footprints, Enable, Source_Type, Source_Name, Fill_Mode, Hatch_Orientation, Hatch_Width, Hatch_Gap, Hatch_Minimum_Hole_Ratio, Smoothing_Effort, Smoothing_Amount, Thieving_Pattern, Thieving_Element_Size, Thieving_Gap, Thieving_Line_Width, Thieving_Stagger, Thieving_Orientation, Remove_Islands, Minimum_Island_Area, Clearance, Minimum_Width, Pad_Connections, Thermal_Relief_Gap, Thermal_Relief_Spoke_Width | Layer -> on_layer(); Name -> group() / prop('Name'); Keep_Out_Tracks -> disallow.kinds track; Keep_Out_Vias -> disallow.kinds via; Keep_Out_Pads -> disallow.kinds pad; Keep_Out_Zone_Fills -> disallow.kinds zone_fill; Keep_Out_Footprints -> disallow.kinds footprint; Clearance -> clearance [item_override] on zones; Minimum_Width -> zone_min_width (value); Pad_Connections -> zone_connection [item_default]; Thermal_Relief_Gap -> thermal_relief_gap [item_default]; Thermal_Relief_Spoke_Width -> thermal_spoke_width [item_default] (all three checked after rules, P9) | pcbnew/zone.cpp:2390† (first registration) |
