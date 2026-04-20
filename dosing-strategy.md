# Aminoglycoside Extended-Interval (EI) Dosing Algorithm

---

## 1. Eligibility for Extended-Interval (EI) Dosing

Patients are eligible for EI dosing unless any exclusion criteria below are present.

---

## 2. Exclusion Criteria for EI Dosing

- Gram-positive synergy indication (e.g., endocarditis)
- Renal insufficiency requiring conventional dosing consideration
- Hemodialysis (HD)
- Continuous renal replacement therapy (CRRT)
- Acute kidney injury with unstable renal function
- Surgical prophylaxis
- Pregnancy or neonatal population
- Nontuberculous mycobacterial (NTM) infection

---

## 3. Dosing Pathways

- **Gram-positive synergy (e.g., endocarditis)** → Conventional dosing
- **Renal insufficiency / AKI** → Conventional dosing
- **Hemodialysis (HD)** → Refer to HD dosing protocol
- **CRRT** → Refer to CRRT dosing protocol
- **Surgical prophylaxis** → Refer to surgical prophylaxis guideline
- **Neonatal population** → Refer to neonatal dosing section
- **NTM infections** → Refer to NTM dosing section

---

## 4. Flowchart

```mermaid
flowchart TD

A[Start: Aminoglycoside Dosing Evaluation] --> B{Eligible for Extended-Interval (EI) Dosing?}

B --> C{Any exclusion criteria present?}

C -- No --> D[Use Extended-Interval Dosing]

C -- Yes --> E{Which exclusion applies?}

E -- Gram-positive synergy --> F[Conventional dosing (e.g., endocarditis)]

E -- Renal insufficiency / AKI --> F

E -- Hemodialysis --> G[Refer to HD dosing protocol]

E -- CRRT --> H[Refer to CRRT dosing protocol]

E -- Surgical prophylaxis --> I[Refer to surgical prophylaxis guidance]

E -- Neonatal --> J[Refer to neonatal dosing section]

E -- NTM infection --> K[Refer to NTM dosing section]
```
