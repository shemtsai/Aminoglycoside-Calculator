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
    Start([Start: Aminoglycoside Evaluation]) --> Criteria

    subgraph Decision_Column [ ]
        direction TB
        Criteria[<b>Exclusion Criteria Present?</b><br/>- Gram-positive synergy<br/>- Renal insufficiency / AKI<br/>- Hemodialysis / CRRT<br/>- Unstable renal function<br/>- Surgical prophylaxis<br/>- Pregnancy / Neonatal<br/>- NTM infection]
        
        %% The Yes path drops vertically
        ExclusionType{Identify Type}
        
        %% The No path also drops vertically
        EI[<b>Extended-Interval Dosing</b>]
    end

    %% Connections
    Criteria -- "Yes" --> ExclusionType
    Criteria -- "No" --> EI

    %% Side Branches from the Identify Type diamond
    ExclusionType -- "All other<br/>(Synergy/AKI)" --> F[<b>Conventional Dosing</b>]
    
    ExclusionType --> G
    ExclusionType --> H
    ExclusionType --> I
    ExclusionType --> J
    ExclusionType --> K

    %% Styling
    style Criteria fill:#fff9c4,stroke:#fbc02d,text-align:left
    style EI fill:#d4edda,stroke:#28a745,stroke-width:2px
    style F fill:#f8d7da,stroke:#dc3545,stroke-width:2px
    style Decision_Column fill:none,stroke:none

    subgraph Specialty_Protocols [Specialty Dosing Protocols]
        direction TB
        G[Hemodialysis]
        H[CRRT]
        I[Surgical Prophylaxis]
        J[Neonatal]
        K[NTM Section]
    end
```
