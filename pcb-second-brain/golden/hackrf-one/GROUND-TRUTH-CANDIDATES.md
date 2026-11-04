# Ground-truth candidates — hackrf-one

Status: **unjudged candidates** produced in W1 (cloud session 2). Each item is a
lead for a revision-diff ground-truth finding (docs/08 § 8.1). Confirmation needs
the kernel's DesignIR diff of the two revisions plus a human look; do not treat
any line below as a confirmed finding.

- Repository: https://github.com/greatscottgadgets/hackrf
- Pinned commit: `468c5d93b66032c892c95bbf9ff15ff1b803e84f` (2026-10-07T10:09:07+02:00)
- Paths analysed: `hardware/hackrf-one`
- Domain: rf
- Commits touching these paths: 234; touching design files: 218; design-file commits whose message matches a fix keyword: 79

## Revision pairs

No revision directories; revisions exist only as git history (pairs = consecutive commits below).

Tags in repository (18): `cft-milestone-1`, `jawbreaker-20120914`, `jawbreaker-20121016`, `2013.06.1`, `v2013.07.1`, `v2014.04.1`, `v2014.08.1`, `v2015.07.1`, `v2015.07.2`, `v2017.02.1`, `v2018.01.1`, `v2021.03.1`, `v2022.09.1`, `v2023.01.1`, `v2024.02.1`, `v2026.01.1`, `v2026.01.2`, `v2026.01.3`

## Revision notes and errata in the repository

Source: `docs/source/list_of_hardware_revisions.rst` at `468c5d93b6` (verbatim)

```text
================================================
Hardware Revisions
================================================
Hardware revisions exist mainly to deal with changes in component availability. Each revision of a product meets the same performance specifications that are measured in the factory.
HackRF Pro
~~~~~~~~~~
The initial production revision of HackRF Pro is r1.2.1.
HackRF One
~~~~~~~~~~
HackRF One r1–r4
^^^^^^^^^^^^^^^^
The first revision of HackRF One shipped by Great Scott Gadgets starting in 2014 was labeled r1. Subsequent manufacturing runs incremented the revision number up to r4 without modification to the hardware design. Manufacturing years: 2014–2020
HackRF One r5
^^^^^^^^^^^^^
This experimental revision has not been manufactured.
HackRF One r6
^^^^^^^^^^^^^
SKY13350 RF switches were replaced by SKY13453. Although the SKY13453 uses simplified control logic, it did not require a firmware modification. Hardware revision detection pin straps were added. Manufacturing year: 2020
HackRF One r7
^^^^^^^^^^^^^
SKY13453 RF switches were reverted to SKY13350. USB VBUS detection resistor values were updated. Manufacturing year: 2021
HackRF One r8
^^^^^^^^^^^^^
SKY13350 RF switches were replaced by SKY13453. Manufacturing years: 2021–2022
HackRF One r9
^^^^^^^^^^^^^
MAX2837 was replaced by MAX2839. Si5351C was replaced by Si5351A with additional clock distribution. A series diode was added to the antenna port power supply. Manufacturing year: 2023
HackRF One r10
^^^^^^^^^^^^^^
This revision is based on r8, reverting most of the changes made in r9. A series diode was added to the antenna port power supply. Manufacturing year: 2024
Hardware Revision Identification
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
HackRF Ones manufactured by Great Scott Gadgets have the revision number printed on the PCB top silkscreen layer near the MAX5864 (U18).
Starting with HackRF One r6, hardware revisions are detected by firmware and reported by ``hackrf_info``.
```

## Issue-tracker leads

GitHub issues matching hardware words (search API, 2026-11-04); titles only, unjudged.

| # | State | Opened | Title |
|---|---|---|---|
| [407](https://github.com/greatscottgadgets/hackrf/issues/407) | open | 2017-09-21 | [Enhancement] Non-streaming transmit and receive in libhackrf? |
| [444](https://github.com/greatscottgadgets/hackrf/issues/444) | open | 2018-01-09 | Hardware documentation improvements |
| [457](https://github.com/greatscottgadgets/hackrf/issues/457) | closed | 2018-02-15 | Install API in Windows |
| [509](https://github.com/greatscottgadgets/hackrf/issues/509) | closed | 2018-08-07 | Suggested LNA/PA replacement model number |
| [686](https://github.com/greatscottgadgets/hackrf/issues/686) | closed | 2020-01-05 | need bin file for U20 W25Q80BV |
| [1061](https://github.com/greatscottgadgets/hackrf/issues/1061) | closed | 2022-03-14 | Fast switching between TX and RX |
| [1392](https://github.com/greatscottgadgets/hackrf/issues/1392) | closed | 2023-12-22 | HackRF One: Replace micro-USB with USB type C |
| [1407](https://github.com/greatscottgadgets/hackrf/issues/1407) | closed | 2024-01-18 | How to set/check gain? |
| [1458](https://github.com/greatscottgadgets/hackrf/issues/1458) | closed | 2024-05-26 | PPTC for LNA protection? |
| [1463](https://github.com/greatscottgadgets/hackrf/issues/1463) | closed | 2024-06-22 | Replacement of TPS62410 Buck converter  |
| [1465](https://github.com/greatscottgadgets/hackrf/issues/1465) | closed | 2024-07-04 | Add USB-c (pcb already modified) |
| [1481](https://github.com/greatscottgadgets/hackrf/issues/1481) | closed | 2024-09-18 | Unwanted Carrier generated at central frequency in TX mode |
| [1494](https://github.com/greatscottgadgets/hackrf/issues/1494) | open | 2024-10-24 | trigger input broken |
| [1502](https://github.com/greatscottgadgets/hackrf/issues/1502) | closed | 2024-11-08 |  check point measurment P15 don't have voltage |
| [1515](https://github.com/greatscottgadgets/hackrf/issues/1515) | closed | 2024-12-13 | USB port not working |
| [1541](https://github.com/greatscottgadgets/hackrf/issues/1541) | closed | 2025-03-28 | How to flash the firmware onto the actual chips? |
| [1560](https://github.com/greatscottgadgets/hackrf/issues/1560) | closed | 2025-06-09 | Lighting Strike |
| [1562](https://github.com/greatscottgadgets/hackrf/issues/1562) | open | 2025-06-21 | Hammond Box 1455J1201 Endplate Template |
| [1573](https://github.com/greatscottgadgets/hackrf/issues/1573) | closed | 2025-08-02 | Hardware upgrade path from HackRF to HackRF Pro |
| [1577](https://github.com/greatscottgadgets/hackrf/issues/1577) | closed | 2025-08-10 | Equivalent or better replacement for the MGA-81563 amplifier LNA chip because its obsolete |
| [1758](https://github.com/greatscottgadgets/hackrf/issues/1758) | open | 2026-05-18 | where are hackrf pro gerbers and pcb files? |
| [1763](https://github.com/greatscottgadgets/hackrf/issues/1763) | closed | 2026-05-25 | hack rf tool build fails |
| [1822](https://github.com/greatscottgadgets/hackrf/issues/1822) | closed | 2026-10-09 | HackRF Pro: hackrf_get_clkin_status() right after hackrf_start_rx() corrupts Si5351 MS0 (I2C read in USB ISR r |

## Commits that look like design fixes

Design-file commits whose subject matches a fix keyword. Categories are keyword
guesses for triage, not judgements. Pair = this commit vs its parent on the same path.

| Commit | Date | Categories | Subject | Design files touched |
|---|---|---|---|---|
| [`d7b5a3537`](https://github.com/greatscottgadgets/hackrf/commit/d7b5a35371194ccb9ecd493bc116b3bc5998daaf) | 2024-02-07 | generic fix | HackRF One: correct fabrication attributes and DNP | baseband.kicad_sch, frontend.kicad_sch, hackrf-one.kicad_pcb, hackrf-one.kicad_sch, mcu.ki |
| [`1d138898a`](https://github.com/greatscottgadgets/hackrf/commit/1d138898ab61bcf8d43fc124b6df15ca3de4c933) | 2024-02-07 | silkscreen | HackRF One: update 0402 silkscreen for readability | hackrf-one.kicad_pcb |
| [`f90930f4c`](https://github.com/greatscottgadgets/hackrf/commit/f90930f4c02777db7fade60ff45436609b51ea3c) | 2024-02-07 | footprint, silkscreen | HackRF One: remove silkscreen from pads | hackrf-one.kicad_pcb |
| [`375ba4e21`](https://github.com/greatscottgadgets/hackrf/commit/375ba4e21ad91dcdbed3dc49526191eabffc4e80) | 2024-02-07 | part swap | HackRF One: clean up symbol value names | baseband.kicad_sch, frontend.kicad_sch, hackrf-one.kicad_pcb, hackrf-one.kicad_pro, hackrf |
| [`44e3b7464`](https://github.com/greatscottgadgets/hackrf/commit/44e3b746473fffb905ab20a2090ee529eae8253b) | 2024-02-07 | generic fix | HackRF One: correct T3 and T4 balun symbol | baseband.kicad_sch, hackrf-one.kicad_pcb |
| [`a4e07bde8`](https://github.com/greatscottgadgets/hackrf/commit/a4e07bde833c1315938a6e53657d6bf19133bfb4) | 2024-02-07 | decoupling, part swap | HackRF One: update 100 nF capacitor part selection | baseband.kicad_sch, frontend.kicad_sch, hackrf-one.kicad_pcb, mcu.kicad_sch |
| [`3f34a678a`](https://github.com/greatscottgadgets/hackrf/commit/3f34a678a5d54657a93de6b7e1b3ab364a823d31) | 2024-02-07 | clearance | HackRF One: add SMA connector keepouts | hackrf-one.kicad_pcb |
| [`ae805c43c`](https://github.com/greatscottgadgets/hackrf/commit/ae805c43c6106a91bf16fa8ad846eeb926f2bddf) | 2023-01-06 | generic fix | HackRF One: correct RXBB net names | baseband.kicad_sch, hackrf-one.kicad_pcb |
| [`9b1ddf8a8`](https://github.com/greatscottgadgets/hackrf/commit/9b1ddf8a8f081e43a27c533c39365a7c4f1d6e9c) | 2023-01-06 | termination | HackRF One: update clock termination passives | clock.kicad_sch, frontend.kicad_sch, hackrf-one.kicad_pcb, hackrf-one.kicad_sch |
| [`220099e3d`](https://github.com/greatscottgadgets/hackrf/commit/220099e3d18e7397f69ffe9a8175dc6ee01de667) | 2023-01-06 | power | HackRF One: update IF balun tuning inductors | baseband.kicad_sch, hackrf-one.kicad_pcb, hackrf-one.kicad_sch |
| [`04f7487a2`](https://github.com/greatscottgadgets/hackrf/commit/04f7487a2f146c10c7b2eec9cdddf0b54d1d276f) | 2023-01-06 | decoupling, part swap | HackRF One: update 100 nF capacitor part selection | baseband.kicad_sch, clock.kicad_sch, frontend.kicad_sch, hackrf-one.kicad_pcb, mcu.kicad_s |
| [`643336b22`](https://github.com/greatscottgadgets/hackrf/commit/643336b224ddaef9be873b53ad8288749e76d871) | 2023-01-06 | clearance | HackRF One: add SMA connector keepouts | hackrf-one.kicad_pcb |
| [`155a4532f`](https://github.com/greatscottgadgets/hackrf/commit/155a4532f0594b9020d65d7433da58f34aba4467) | 2023-01-06 | part swap | HackRF One: swap pins for platform detection | hackrf-one.kicad_pcb, mcu.kicad_sch |
| [`895be8964`](https://github.com/greatscottgadgets/hackrf/commit/895be8964f560b3111cd57ae6c081025ea84bcc1) | 2023-01-06 | decoupling | HackRF One: add diode bypass jumper | frontend.kicad_sch, hackrf-one.kicad_pcb, hackrf-one.kicad_sch |
| [`793c366dc`](https://github.com/greatscottgadgets/hackrf/commit/793c366dcd86eaca8cdaaf7cdc05920f4df5cdfa) | 2023-01-06 | generic fix | HackRF One: correct or ignore ERC errors | baseband.kicad_sch, hackrf-one.kicad_pro, hackrf-one.kicad_sch, mcu.kicad_sch |
| [`2c456a5a1`](https://github.com/greatscottgadgets/hackrf/commit/2c456a5a16f271ffd3b88d8b21922afa3a1d7904) | 2023-01-06 | routing | HackRF One: update trace widths for new stack-up | baseband.kicad_sch, frontend.kicad_sch, hackrf-one.kicad_pcb, hackrf-one.kicad_pro, hackrf |
| [`fb48e629c`](https://github.com/greatscottgadgets/hackrf/commit/fb48e629c18b7204ae22e94a74cc059da1985345) | 2021-09-27 | part swap | HackRF One: update VBUS sense resistor values | hackrf-one.kicad_pcb, mcu.sch |
| [`7c264b90c`](https://github.com/greatscottgadgets/hackrf/commit/7c264b90cbaf5054a1e469b42d08da2477f648a7) | 2020-10-13 | silkscreen | HackRF One: update silkscreen | hackrf-one.kicad_pcb |
| [`585e69a5d`](https://github.com/greatscottgadgets/hackrf/commit/585e69a5d06ffdf2077f33d62991dfa6523dcd12) | 2020-10-13 | part swap | HackRF One: replace obsolete SKY13350-385LF RF switch with SKY13453-385LF | baseband.sch, frontend.sch, hackrf-one-cache.lib, hackrf-one.kicad_pcb, mcu.sch |
| [`5e9cad663`](https://github.com/greatscottgadgets/hackrf/commit/5e9cad66363467fae93aec29679c041c03905047) | 2018-04-30 | part swap | HackRF One: replaced U20 with Winbond W25Q80DVSSIG | baseband.sch, frontend.sch, hackrf-one.kicad_pcb, hackrf-one.pro, hackrf-one.sch, mcu.sch |
| [`b4a452cd3`](https://github.com/greatscottgadgets/hackrf/commit/b4a452cd3a6e77c3c6a7ac324566c82905247e5f) | 2018-04-14 | decoupling, part swap | replaced AVX/Kyocera 2500BL14M100T (NRND) with Abracon ABM8-25.000MHZ-10-D1G-T (X1), added load capacitors (C1 | baseband.sch |
| [`26360e18f`](https://github.com/greatscottgadgets/hackrf/commit/26360e18fdca595688d5740a5a1be2c43134e043) | 2018-04-14 | part swap | listed Murata LDB182G4510G-120 as an alternative to Johanson 2500BL14M100T (T3 and T4) | baseband.sch |
| [`86f724a34`](https://github.com/greatscottgadgets/hackrf/commit/86f724a342fb2e0d8423674083b9a046a82e6d3d) | 2018-04-14 | part swap | replaced Murata LXES1TBCC2-004 (End of Life) with Vishay VBUS54CV-HSF-G4-08 (U15) | mcu.sch |
| [`370ee2dce`](https://github.com/greatscottgadgets/hackrf/commit/370ee2dce5c5d26753ca7f725c2d1d5dbf832d1e) | 2018-02-28 | part swap | HackRF One: replaced U20 with Winbond W25Q80JVSSIQ | mcu.sch |
| [`576e25b1d`](https://github.com/greatscottgadgets/hackrf/commit/576e25b1d4f4d0d544e674f9702e70fa4c48fea3) | 2014-05-19 | silkscreen | updated date on silk screen | hackrf-one.brd |
| [`43f96036e`](https://github.com/greatscottgadgets/hackrf/commit/43f96036e4efb51e1fbaecb5115f640c36387303) | 2014-05-19 | footprint | increased SKY13317 pad size | hackrf-one.brd |
| [`cff03762e`](https://github.com/greatscottgadgets/hackrf/commit/cff03762e7c8d97497e1607e585fadbafbdaa49c) | 2014-03-06 | generic fix | fixed USB port position on end plate | hackrf-one.brd |
| [`e78a3f9ed`](https://github.com/greatscottgadgets/hackrf/commit/e78a3f9edc36939a26fc9afbb84e3d595eb1e4c2) | 2014-02-13 | part swap, power | H1 replaced 1 uH inductors with 10 uH | frontend.sch, hackrf-one-cache.lib, hackrf-one.brd |
| [`5d7770f17`](https://github.com/greatscottgadgets/hackrf/commit/5d7770f178afed219d3d84d82e4660a5016ba800) | 2014-02-13 | decoupling, part swap | H1 replaced 47 pF IF DC blocking caps with 22 pF | frontend.sch, hackrf-one-cache.lib, hackrf-one.brd |
| [`aaf5258cb`](https://github.com/greatscottgadgets/hackrf/commit/aaf5258cb00a0b13bb89b747cfc563cffd34bc93) | 2014-02-13 | decoupling, part swap | H1 replaced 10 nF DC blocking caps with 100 nF | frontend.sch, hackrf-one.brd |
| [`580edaa53`](https://github.com/greatscottgadgets/hackrf/commit/580edaa533b08d37b092a76bfa7e4fe73b1ac05b) | 2014-02-13 | generic fix | H1 fixed power supply switch by removing Q6 | baseband.sch, frontend.sch, hackrf-one-cache.lib, hackrf-one.brd, hackrf-one.sch, mcu.sch |
| [`db4221cac`](https://github.com/greatscottgadgets/hackrf/commit/db4221caca73077ac5a3f1eead13a9109ac117a1) | 2014-02-11 | routing | H1 made power supply vias bigger | hackrf-one.brd |
| [`cede949d5`](https://github.com/greatscottgadgets/hackrf/commit/cede949d5b53593654428df52af49ca37b5eb5b6) | 2014-02-04 | part swap | restored 0 ohm resistors on USB signals, could be replaced by DLP11S common mode choke | mcu.sch |
| [`45115b6ab`](https://github.com/greatscottgadgets/hackrf/commit/45115b6abc93e0fb53fcb0961d220932bfa1061e) | 2014-01-20 | generic fix | fixed INSERT attributes | hackrf-one.brd |
| [`f6c153e73`](https://github.com/greatscottgadgets/hackrf/commit/f6c153e73d0331434b8c47ea9e0df52f96f71405) | 2014-01-11 | generic fix | switch to fixed SMA module | hackrf-one.brd |
| [`ceb49ae0a`](https://github.com/greatscottgadgets/hackrf/commit/ceb49ae0ae4175e8ee2f56ce395826ee3ec0c79b) | 2014-01-10 | part swap | bom update for new SMA connectors | baseband.sch, frontend.sch, hackrf-one-cache.lib, mcu.sch |
| [`a32724c11`](https://github.com/greatscottgadgets/hackrf/commit/a32724c1193befe27afbb21dd8c614aa62c3706e) | 2013-12-12 | routing | removed a couple GND pins to shorten P20, improving laser cut enclosure ability | hackrf-one-cache.lib, hackrf-one.brd, mcu.sch |
| [`72662dbbc`](https://github.com/greatscottgadgets/hackrf/commit/72662dbbc4a4c81d05e99e44fca436a5ef071ee8) | 2013-12-12 | decoupling, footprint | load cap footprints for main crystal just in case | baseband.sch, frontend.sch, hackrf-one-cache.lib, hackrf-one.brd, mcu.sch |
| [`ba7f3c86e`](https://github.com/greatscottgadgets/hackrf/commit/ba7f3c86ea7b647266d2ea9d1a5d491ed3f437bb) | 2013-12-12 | pull-up/down, part swap | changed some 4k7 pull-up or pull-down resistors to 10k to reduce BOM | baseband.sch, frontend.sch, hackrf-one-cache.lib, hackrf-one.brd |
| [`a07762bb3`](https://github.com/greatscottgadgets/hackrf/commit/a07762bb3e41cfab3ebeef6fae754af684e3878e) | 2013-12-12 | generic fix | R16 description corrected | baseband.sch, frontend.sch, hackrf-one-cache.lib, hackrf-one.sch, mcu.sch |
| [`669c02712`](https://github.com/greatscottgadgets/hackrf/commit/669c027122c8edc39223156d18144db19d40eb9a) | 2013-12-12 | silkscreen | PCB layer markings in copper | hackrf-one.brd |
| [`0eb2e9435`](https://github.com/greatscottgadgets/hackrf/commit/0eb2e9435d91398314aee4cc93b056fdd4c4c22c) | 2013-11-29 | decoupling, part swap | replaced obsolete capacitor selections | baseband.sch, frontend.sch, hackrf-one.brd, mcu.sch |
| [`91b42ed9c`](https://github.com/greatscottgadgets/hackrf/commit/91b42ed9c6597d58f56c02c3ec383dd45b6184e9) | 2013-11-27 | silkscreen | date on silkscreen | hackrf-one.brd |
| [`dd7b3a2da`](https://github.com/greatscottgadgets/hackrf/commit/dd7b3a2dade917f6e29188c1bb3a3d9d4cc38afa) | 2013-11-27 | decoupling | SSP1_MISO decoupling cap just in case | baseband.sch, hackrf-one-cache.lib, hackrf-one.brd, mcu.sch |
| [`c9f0efa28`](https://github.com/greatscottgadgets/hackrf/commit/c9f0efa288156024db3740c3368581cc37c6ad8f) | 2013-11-27 | routing, generic fix | trace bug introduced previous commit | hackrf-one.brd |
| [`f00747f26`](https://github.com/greatscottgadgets/hackrf/commit/f00747f260695ecc30dd4345b6ecd670ba99aefc) | 2013-11-25 | silkscreen | relabeled LED2 and LED3 to RX and TX | hackrf-one-cache.lib, hackrf-one.brd, mcu.sch |
| [`eab5418a5`](https://github.com/greatscottgadgets/hackrf/commit/eab5418a51eb1ddc42d258e8a597c2d54ab7819d) | 2013-11-25 | routing | optional passives between USB shield and GND | baseband.sch, frontend.sch, hackrf-one-cache.lib, hackrf-one.brd, mcu.sch |
| [`693935cd0`](https://github.com/greatscottgadgets/hackrf/commit/693935cd0d72723ca0f78febd17ec4f8200e9e4a) | 2013-11-25 | routing | clippable GND test point | baseband.sch, frontend.sch, hackrf-one-cache.lib, hackrf-one.brd, hackrf-one.sch, mcu.sch |
| [`cd96c356e`](https://github.com/greatscottgadgets/hackrf/commit/cd96c356e2aad81c4400ff6892aaec24bcfaaacd) | 2013-11-22 | footprint, protection, routing | lengthened USB TVS pads | hackrf-one.brd |
| [`575a8394e`](https://github.com/greatscottgadgets/hackrf/commit/575a8394ed35facc1001faa1ee16b6da4e05068b) | 2013-11-22 | generic fix | fixed GP_CLKIN passives | baseband.sch, frontend.sch, hackrf-one-cache.lib, hackrf-one.brd, hackrf-one.sch, mcu.sch |
| [`15a51cab5`](https://github.com/greatscottgadgets/hackrf/commit/15a51cab55ddc7b33862ff841da5405681747dcf) | 2013-11-22 | routing | grounded pushbutton mechanical support | baseband.sch, frontend.sch, hackrf-one-cache.lib, hackrf-one.brd, mcu.sch |
| [`04898a782`](https://github.com/greatscottgadgets/hackrf/commit/04898a78205d383fbeda8d22783298c54b65e65a) | 2013-11-22 | routing | exposed GND through mounting holes | baseband.sch, frontend.sch, hackrf-one-cache.lib, hackrf-one.brd, mcu.sch |
| [`0002351b2`](https://github.com/greatscottgadgets/hackrf/commit/0002351b21a2f9f02de4a326b64b6e4a1c6c8f5d) | 2013-11-21 | decoupling, part swap | changed some capacitor values to improve low frequency performance | baseband.sch, frontend.sch, hackrf-one-cache.lib, hackrf-one.brd, hackrf-one.sch, mcu.sch |
| [`95ffc704a`](https://github.com/greatscottgadgets/hackrf/commit/95ffc704a1cd2543db943395cb78881ebe65fde1) | 2013-11-16 | generic fix | P28 and P29 reworked, exposed unused SGPIO signals, moved some CPLD JTAG signals to P28 | baseband.sch, frontend.sch, hackrf-one-cache.lib, hackrf-one.brd, hackrf-one.sch, mcu.sch |
| [`06f345239`](https://github.com/greatscottgadgets/hackrf/commit/06f345239ba9a9f5cf638593a4f7422c470f27c2) | 2013-11-11 | silkscreen | silkscreen tweaks | hackrf-one.brd |
| [`fecc7346b`](https://github.com/greatscottgadgets/hackrf/commit/fecc7346b3c7a5b4d28334207b400578d68f3e53) | 2013-11-11 | routing | GND test points | baseband.sch, frontend.sch, hackrf-one-cache.lib, hackrf-one.brd, mcu.sch |
| [`a8c2c0b6d`](https://github.com/greatscottgadgets/hackrf/commit/a8c2c0b6d17b6aaf4f0b51d3fb9ce3cf33413ed0) | 2013-11-11 | decoupling | more decoupling caps | baseband.sch, frontend.sch, hackrf-one-cache.lib, hackrf-one.brd, mcu.sch |
| [`26104e673`](https://github.com/greatscottgadgets/hackrf/commit/26104e673535554fa57f36e992c76adfe3c6eb81) | 2013-11-11 | routing | nudged some traces in the RF section | hackrf-one.brd |
| [`c5533b3c9`](https://github.com/greatscottgadgets/hackrf/commit/c5533b3c96b5419ca03ac8c3454abac35c719025) | 2013-11-11 | routing, generic fix | reworked zones so LED signals do not cross power planes | baseband.sch, frontend.sch, hackrf-one-cache.lib, hackrf-one.brd, mcu.sch |
| [`f576fc27f`](https://github.com/greatscottgadgets/hackrf/commit/f576fc27f07d39c799b82fe2e2ac5ec530f4fa98) | 2013-11-11 | routing | rerouted 1V8 | hackrf-one.brd |
| [`5beef42bc`](https://github.com/greatscottgadgets/hackrf/commit/5beef42bc59bc57ad594742a53b3491f156f9d4b) | 2013-11-10 | routing | track clean-up | hackrf-one.brd |
| [`198bc7afe`](https://github.com/greatscottgadgets/hackrf/commit/198bc7afeb1611240572b68a3082fe702f8b1b80) | 2013-11-10 | protection | TVS diodes for CLKIN and CLKOUT | baseband.sch, frontend.sch, hackrf-one-cache.lib, hackrf-one.brd, mcu.sch |
| [`61176e537`](https://github.com/greatscottgadgets/hackrf/commit/61176e53718cd73ba51a69aa51aaad0093f3c752) | 2013-11-10 | part swap, protection | selected a new part for USB ESD protection | hackrf-one-cache.lib, hackrf-one.brd, mcu.sch |
| [`02781f6bc`](https://github.com/greatscottgadgets/hackrf/commit/02781f6bc5e68b6d35bfe3419dd97c67e8d1f681) | 2013-11-09 | footprint | extra footprints for clock signal passives just in case | baseband.sch, frontend.sch, hackrf-one-cache.lib, hackrf-one.brd, mcu.sch |
| [`f139288f6`](https://github.com/greatscottgadgets/hackrf/commit/f139288f6206ed61d02f0e6095bb9f0ed02dd297) | 2013-11-09 | generic fix | fixed SMA connector selection | baseband.sch, frontend.sch, hackrf-one.brd |
| [`ee5e6f1df`](https://github.com/greatscottgadgets/hackrf/commit/ee5e6f1dfa6e68874fba73830092419d8759d90e) | 2013-11-09 | footprint | paste layer cleanup | hackrf-one.brd |
| [`08655a1fa`](https://github.com/greatscottgadgets/hackrf/commit/08655a1fa1f5c4de8ddf1269b7ab636af3bc7d9a) | 2013-11-09 | clearance, routing, generic fix | new style PCB trace jumper to minimize DRC errors | baseband.sch, hackrf-one-cache.lib, hackrf-one.brd, mcu.sch |
| [`2a04fa4d1`](https://github.com/greatscottgadgets/hackrf/commit/2a04fa4d10de675f8cfac472c0b38faa44a5a724) | 2013-11-09 | clearance | more DRC cleanup | hackrf-one.brd |
| [`7c03751e2`](https://github.com/greatscottgadgets/hackrf/commit/7c03751e21903d04fb0184cbb530d5eac1ccca27) | 2013-11-09 | clearance | DRC cleanup | hackrf-one.brd |
| [`bef5835d5`](https://github.com/greatscottgadgets/hackrf/commit/bef5835d54aacc740cb80bf727c4639b2cf91759) | 2013-11-01 | power | USB, regulator layout | hackrf-one-cache.lib, hackrf-one.brd, mcu.sch |
| [`7bfd46597`](https://github.com/greatscottgadgets/hackrf/commit/7bfd46597b7054d49544545cc89ba1b430ac5159) | 2013-10-25 | pull-up/down | USB0_ID pull-up | baseband.sch, frontend.sch, hackrf-one-cache.lib, mcu.sch |
| [`0d59261ae`](https://github.com/greatscottgadgets/hackrf/commit/0d59261ae32eeceb0720f4552b23906cb56e6519) | 2013-10-25 | part swap | removed RF switch logic ICs, replaced with direct GPIO | baseband.sch, frontend.sch, hackrf-one-cache.lib, hackrf-one.sch, mcu.sch |
| [`d5fbf3554`](https://github.com/greatscottgadgets/hackrf/commit/d5fbf355454cf2e8fed9f79cb2b92cf27f53e528) | 2013-10-24 | generic fix | fixed all clock generator output supplies to VCC | baseband.sch, hackrf-one-cache.lib |
| [`cbd3295cb`](https://github.com/greatscottgadgets/hackrf/commit/cbd3295cb35ad5ea8d94f5240cf7985024589f06) | 2013-10-24 | generic fix | fixed clock input and output to 3.3V CMOS | baseband.sch, frontend.sch, hackrf-one-cache.lib, hackrf-one.sch, mcu.sch |
| [`40406e52a`](https://github.com/greatscottgadgets/hackrf/commit/40406e52ade66c27dbd318773877a00db97d0d12) | 2013-09-29 | part swap | removed more DNP parts | hackrf-one-cache.lib, mcu.sch |
| [`cda1462ff`](https://github.com/greatscottgadgets/hackrf/commit/cda1462fff831d65fbe4f72fc5c30c3ebb532291) | 2013-09-29 | power, generic fix | regulator feedback fix | hackrf-one-cache.lib, mcu.sch |
| [`2626e1814`](https://github.com/greatscottgadgets/hackrf/commit/2626e1814d5f16e4350495203124ee6ae54297c9) | 2013-09-29 | generic fix | fixed MAX5864 OVDD to VCC | baseband.sch, hackrf-one-cache.lib |
| [`a9026c521`](https://github.com/greatscottgadgets/hackrf/commit/a9026c521f51255eb11e706fe776b2f4780baf80) | 2013-09-29 | generic fix | fixed VCCIO1 to VCC | hackrf-one-cache.lib, mcu.sch |
| [`dae882035`](https://github.com/greatscottgadgets/hackrf/commit/dae88203572fd0ff7785fcf88c2e906da88c8a64) | 2013-09-29 | part swap | removed some DNP parts and 0 ohm resistors | hackrf-one-cache.lib, mcu.sch |
