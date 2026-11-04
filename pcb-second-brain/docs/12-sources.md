# 12 — Sources consulted (November 2026)

These were checked while writing this plan; re-verify at implementation
time since tooling moves fast. Add a "Corrections" section as reality
diverges from the plan.

## KiCad 10, IPC API, Python bindings
- kicad-mcp-kipy on PyPI (describes KiCad 10 IPC API over NNG + `kipy`): https://pypi.org/project/kicad-mcp-kipy/
- KiCad bridge capability notes (schematic IPC in KiCad 11): https://wiki.adom.inc/adom/kicad-bridge/files/capabilities.py
- KiCad Buddy (notes on KiCad 10 IPC schematic limits): https://pypi.org/project/kicad-buddy/
- KiCad forum thread on IPC API and KiCad 10 connection issues: https://forum.kicad.info/t/kicad-ipc-api-attached-to-kicad-9-but-not-to-kicad-10/68030

## Schematic-as-code
- atopile repository (MIT): https://github.com/atopile/atopile
- atopile README: https://cdn.jsdelivr.net/gh/atopile/atopile@main/README.md
- atoplace (third-party layout assistant on KiCad IPC): https://github.com/glovepost/atoplace

## KiCad MCP servers (community; evaluate, do not adopt blindly)
- KiCad forum: "I built an MCP server that lets AI assistants help with KiCad PCB design": https://forum.kicad.info/t/i-built-an-mcp-server-that-lets-ai-assistants-help-with-kicad-pcb-design/66451/8
- lamaalrajih/kicad-mcp (analysis-focused): https://playbooks.com/mcp/lamaalrajih-kicad
- ProductOfAmerica/mcp-server-kicad (100+ tools, early stage): https://enterprisedna.co/directories/mcp/productofamerica-mcp-server-kicad
- oaslananka/kicad-mcp-pro: https://mcp.so/ja/servers/kicad-mcp-pro
- Directory listing: https://lobehub.com/mcp?q=pcb

## Component data APIs
- Nexar API overview: https://nexar.com/api
- Nexar/Octopart data consistency note: https://support.nexar.com/support/solutions/articles/101000541169-differences-between-octopart-com-and-nexar-api-data
- Altium-authored explainer on Nexar API: https://resources.altium.com/p/nexar-api-explained
- Third-party review of Octopart/Nexar tiers (unofficial): https://www.ki-syndikat.de/tools/octopart/

## Routing engines
- Quilter 2026 autorouting review (vendor): https://quilter.ai/blog/pcb-autorouting-in-2026-a-review-of-traditional-tools-vs-quilters-ai-approach
- DeepPCB 2026 comparison (vendor): https://deeppcb.ai/what-is-the-best-ai-pcb-autorouter-in-2026/
- Independent-ish write-up comparing routers: https://dev.to/anypcba_official/i-let-ai-design-a-pcb-heres-how-it-went-5e41
- JLCPCB blog on AI tools by workflow stage: https://jlcpcb.com/blog/ai-tools-for-pcb-design

## Corrections
(none yet)
