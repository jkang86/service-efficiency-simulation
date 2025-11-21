# 📊 Service Efficiency Simulation  
*A Monte Carlo & Taguchi Loss Function Approach*

This repository contains a full simulation study modeling customer wait times, cashier idle times, and financial loss using **Monte Carlo simulation** and the **Taguchi Loss Function**.  
All analysis is performed in **Julia**, using Poisson-distributed arrivals and Normally-distributed service times.

---

## 🚀 Overview

Real-world service environments (retail, banking, customer service centers, etc.) must balance:

- **Customer satisfaction** (short wait times)  
- **Operational efficiency** (minimal cashier idle time)

This project models that trade-off using:

- 🧮 **Probability distributions**  
- 🔁 **Monte Carlo simulation**  
- 💸 **Taguchi Loss Function**  
- 📈 **Visual analysis** (histograms + loss curves)

---

## 📐 Methodology

### 1. Customer Arrivals  
Modeled with a **Poisson distribution**  
\[
\lambda = 10
\]
per unit time  
(PDF Page 2 — Poisson justification) :contentReference[oaicite:2]{index=2}

### 2. Service Times  
Modeled with a **Normal distribution**  
\[
\mu = 7,\quad \sigma = 2
\]
Service times forced positive using `abs()`.

### 3. Monte Carlo Simulation  
- **1,000 trials**  
- Tracks:
  - Customer **wait times**
  - Cashier **idle times**
  - Using queue logic (single-server model)

### 4. Taguchi Loss Function  
Used to measure financial loss from deviations from optimal service time:

\[
L(y) = k(y - m)^2
\]

Based on PDF values (Page 4):  
- Optimal time: **m = 2 minutes**  
- Loss coefficients:
  - **k = 10** for wait time  
  - **k = 5** for idle time  

---

## 🧠 Key Results (from PDF)

From simulation output (PDF Page 3): :contentReference[oaicite:3]{index=3}

### **Wait Time**
- Mean: **0.4234 min**
- SD: **1.0868 min**

### **Idle Time**
- Mean: **3.4149 min**
- SD: **3.1243 min**

### **Taguchi Loss Values**
(PDF Page 4–5)  
- **Wait Time Loss:** 24.79  
- **Idle Time Loss:** 10.02  

Interpretation (PDF Page 6):  
- Long waits → high customer dissatisfaction  
- High idle times → financial inefficiency  
- Taguchi curve shows losses grow **quadratically** with deviation from target

---

## 🧩 Julia Code (Core Implementation)

### 🔹 Monte Carlo Simulation  
(From PDF Page 3)  
```julia
using Random, Distributions, Plots

Random.seed!(42)

lambda_arrival = 10
mu_service = 7
sigma_service = 2
num_trials = 1000

arrival_times = cumsum(rand(Poisson(lambda_arrival), num_trials))
service_times = abs.(rand(Normal(mu_service, sigma_service), num_trials))

wait_times = zeros(num_trials)
idle_times = zeros(num_trials)

server_busy_until = 0.0

for i in 1:num_trials
    if arrival_times[i] < server_busy_until
        wait_times[i] = server_busy_until - arrival_times[i]
    else
        idle_times[i] = arrival_times[i] - server_busy_until
    end
    server_busy_until = arrival_times[i] + service_times[i]
end
