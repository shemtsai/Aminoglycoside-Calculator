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
    Start([Start: Aminoglycoside Evaluation]) --> Criteria{Exclusion Criteria<br/>Present?}
    
    %% Main Path
    Criteria -- No --> EI[<b>Extended-Interval Dosing</b>]
    style EI fill:#d4edda,stroke:#28a745

    %% Exclusion Path
    Criteria -- Yes --> ExclusionType{Exclusion Type}

    %% Grouping Outcomes
    subgraph Specialty_Protocols [Specialty Dosing Protocols]
        direction TB
        G[Hemodialysis Protocol]
        H[CRRT Protocol]
        I[Surgical Prophylaxis]
        J[Neonatal Section]
        K[NTM Section]
    end

    ExclusionType -->|Synergy / AKI| F[<b>Conventional Dosing</b>]
    style F fill:#f8d7da,stroke:#dc3545

    ExclusionType -->|HD| G
    ExclusionType -->|CRRT| H
    ExclusionType -->|Surgery| I
    ExclusionType -->|Neonatal| J
    ExclusionType -->|NTM| K
```
