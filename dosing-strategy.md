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
    Start([Start: Aminoglycoside Evaluation]) --> Criteria{Exclusion Criteria<br/>Present?}
    
    %% The "No" path drops vertically
    subgraph Main_Path [ ]
        direction TB
        Criteria
        EI[<b>Extended-Interval Dosing</b>]
    end
    Criteria -- No --> EI
    style EI fill:#d4edda,stroke:#28a745,stroke-width:2px
    style Main_Path fill:none,stroke:none

    %% The "Yes" path leads to the List, then the List leads to Identify Type
    Criteria -- Yes --> List[Exclusion Criteria:<br/>- Gram-positive synergy<br/>- Renal insufficiency / AKI<br/>- Hemodialysis / CRRT<br/>- Unstable renal function<br/>- Surgical prophylaxis<br/>- Pregnancy / Neonatal<br/>- NTM infection]
    style List fill:#fff9c4,stroke:#fbc02d,text-align:left

    %% Identify Type now branches from the List box
    List --> ExclusionType{Identify Type}

    ExclusionType --> F[<b>Conventional Dosing</b>]
    style F fill:#f8d7da,stroke:#dc3545,stroke-width:2px

    %% Specialty protocols on the far right
    subgraph Specialty_Protocols [Specialty Dosing Protocols]
        direction TB
        G[Hemodialysis]
        H[CRRT]
        I[Surgical Prophylaxis]
        J[Neonatal]
        K[NTM Section]
    end

    ExclusionType --> G
    ExclusionType --> H
    ExclusionType --> I
    ExclusionType --> J
    ExclusionType --> K
```
