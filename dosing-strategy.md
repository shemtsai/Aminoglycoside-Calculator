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
    %% Step 1: Using spacers to force a wide rectangle
    Start["&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>Step 1: EI Eligibility Checklist</b>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<br/>- Gram-positive synergy<br/>- Renal insufficiency / AKI / CrCl 20<br/>- Hemodialysis / CRRT<br/>- Surgical prophylaxis<br/>- Pregnancy / Neonatal<br/>- NTM infection"] --> B{Any Exclusions<br/>Present?}
    
    style Start fill:#fff9c4,stroke:#fbc02d,text-align:left

    %% The 'No' Path
    B -- "No" --> EI["&nbsp;&nbsp;&nbsp;<b>Extended-Interval Dosing</b>&nbsp;&nbsp;&nbsp;"]
    style EI fill:#d4edda,stroke:#28a745,stroke-width:2px

    %% The 'Yes' Path
    B -- "Yes" --> C{"Is it Synergy, AKI,<br/>or CrCl 20?"}
    
    C -- "Yes" --> Conv["&nbsp;&nbsp;&nbsp;<b>Conventional Dosing</b>&nbsp;&nbsp;&nbsp;"]
    style Conv fill:#f8d7da,stroke:#dc3545,stroke-width:2px

    %% Step 3: Wide referral box
    C -- "No" --> D["&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>Refer to Specific Dosing Section</b>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<br/>or Contact ID Pharmacy<br/><br/>- Hemodialysis<br/>- CRRT<br/>- Surgical Prophylaxis<br/>- Neonatal Population<br/>- NTM Infections"]
    
    style D fill:#e1f5fe,stroke:#01579b,text-align:left
```
