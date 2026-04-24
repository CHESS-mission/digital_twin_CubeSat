# CubeSat EST - Spacecraft Parameters

Explains the key choices made in `spacecraft_template.json`.

**Sources:** `CP0_Power_Budget_2026-03-22.xlsx` (sheets: ADCS-New, Transceiver-New, GNSS-New, EPS-New, Total-New) and `CP0_Energy_Budget_2026-03-22.xlsx`.

**Margin policy:** ADCS, GNSS, EPS: +25% | Transceiver: +5%

---

## EPS - Battery

- `total_energy`: **50 Wh** - updated battery capacity from CP0 Energy Budget reference
- `min_battery`: **20% of total** - same ratio as previous design
- `max_battery`: 100% of total - simulation charges to full
- `init_battery_level`: 100% - simulation starts with full battery
- `measure_threshold`: 40% of total - enough energy to run a full measurement campaign
- `com_threshold`: 45% of total - slightly above measurement threshold as UHF com consumes more
- `xb_threshold`: 60% of total - X-band consumes ~30 W (highest mode), requires a well-charged battery

## EPS Board

- `consumption`: 0.3034 W (EPS-New sheet avg) × 1.25 = **0.379 W**, constant across all modes

## Solar Panels

- 28 cells (4 panels × 7 cells), unchanged from previous design
- `efficiency`: 0.25 (BOL) - unchanged 

---

## Telecom

All values from `Transceiver-New` sheet × 1.05 margin.

---

## ADCS

All values from `ADCS-New` sheet (CubeSpace) × 1.25 margin. Mode-dependent:

---

## OBC

- `consumption`: 1.7 W - unchanged from previous design 

---

## Payload

- `consumption_gnss`: GNSS-New Continuous Tracking × 1.25, active in all modes except SAFE
- `consumption_measurment`: represents **Novoviz instrument** active consumption (mode 1 only), from Total-New sheet
