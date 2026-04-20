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
    %% Step 1: Initial Screen
    A[<b>Step 1: EI Eligibility</b><br/>Check Exclusion List:<br/>- Synergy/AKI/Low CrCl<br/>- HD/CRRT<br/>- Pregnancy/Neonatal<br/>- Surgery/NTM] --> B{Any Exclusions<br/>Present?}
    style A fill:#fff9c4,stroke:#fbc02d,text-align:left

    %% Path to EI
    B -- No --> EI[<b>Extended-Interval Dosing</b>]
    style EI fill:#d4edda,stroke:#28a745,stroke-width:2px

    %% Step 2: The Major Exclusions
    B -- Yes --> C{Is it Synergy, AKI,<br/>or CrCl < 20?}
    
    C -- Yes --> Conv[<b>Conventional Dosing</b>]
    style Conv fill:#f8d7da,stroke:#dc3545,stroke-width:2px

    %% Step 3: The Specialty Referrals
    C -- No --> D[<b>Refer to Specific Dosing Section</b><br/><i>or Contact ID Pharmacy</i>]
    style D fill:#e1f5fe,stroke:#01579b

    %% Branching out to specialty nodes for clarity
    direction LR
    D --> G[Hemodialysis]
    D --> H[CRRT]
    D --> I[Surgical Prophylaxis]
    D --> J[Neonatal]
    D --> K[NTM Section]
```
