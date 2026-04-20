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
flowchart LR
    %% Main Entry Point
    Criteria{<b>Are any Exclusion<br/>Criteria Present?</b>}
    
    %% The "Sticky Note" List
    Criteria --- List[<b>Exclusion Criteria:</b><br/>• Gram-positive synergy<br/>• Renal insufficiency / AKI<br/>• Hemodialysis / CRRT<br/>• Unstable renal function<br/>• Surgical prophylaxis<br/>• Pregnancy / Neonatal<br/>• NTM infection]
    style List fill:#fff9c4,stroke:#fbc02d,text-align:left

    %% The "No" Path (Directly to EI)
    Criteria -- No --> EI[<b>Extended-Interval Dosing</b>]
    style EI fill:#d4edda,stroke:#28a745,stroke-width:2px

    %% The "Yes" Path (Directly to specific categories)
    Criteria -- Yes --> F[<b>Conventional Dosing</b><br/><i>(Synergy / AKI)</i>]
    style F fill:#f8d7da,stroke:#dc3545,stroke-width:2px

    %% Direct links to specialty protocols
    subgraph Specialty_Protocols [Specialty Dosing Protocols]
        direction TB
        G[Hemodialysis]
        H[CRRT]
        I[Surgical Prophylaxis]
        J[Neonatal]
        K[NTM Section]
    end

    Criteria -- Yes --> G
    Criteria -- Yes --> H
    Criteria -- Yes --> I
    Criteria -- Yes --> J
    Criteria -- Yes --> K
```
