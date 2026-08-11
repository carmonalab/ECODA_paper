# Source: https://doc.eresearch.unige.ch/hpc/PI_Resources
# Snapshot: 2026-08-11
# Crawled: 2026-08-11T14:32:25Z

---

## Primary Investigator (PI) resources information

A High Performance Computing (HPC) cluster is a shared scientific infrastructure that provides a large amount of computing resources for research projects.

At the University of Geneva, the HPC infrastructure is operated as the **Baobab HPC Service** and consists of three HPC clusters: **Baobab**, **Yggdrasil**, and **Bamboo**.

As Baobab was the first HPC cluster deployed at the University and remains the most well-known among our users, we chose to use its name to designate the entire HPC service.

This page is intended for Principal Investigators (PIs) and explains how **Billing Units** (BU), our unified measure of computing resource usage, are allocated to research groups, how consumption is monitored, the costs associated with HPC usage, and the options available for **purchasing** or **renting** private compute nodes.
---

### Cost Model

The HPC service is free in specific situations:

- **Educational courses** can use the platform free of charge.
- Each Principal Investigator (PI) receives an annual allocation of free **Billing Units (BU)**.

Additional usage can be covered through:

- **[Pay-per-use billing](https://doc.eresearch.unige.ch/pi_resources#price_per_billing_unit)** according to the SNSF pricing table.
- **[Purchasing or renting compute nodes](https://doc.eresearch.unige.ch/pi_resources#purchasing_or_renting_private_compute_nodes)**, which provides additional annual Billing Units and higher scheduling priority.

#### Cost Model at a Glance

The HPC accounting model is based on a single resource unit called the **Billing Unit (BU)**.

- Every PI receives **100,000 BU per year** free of charge.
- Research groups that own or rent compute nodes receive additional annual BU allocations.
- CPU, memory, and GPU usage are converted into BU using standardized conversion factors.
- BU can be consumed on any of our HPC clusters: **Baobab**, **Yggdrasil**, and **Bamboo**.
- Once all available BU have been consumed, additional usage is billed according to the SNSF pricing table.

#### Annual Billing Units

Research groups may receive Billing Units from two sources:

- **100,000 Billing Units annual** free allocation provided to every PI.
- The hardware owned or rented by the research group.

You can check your available BU using:

[Resources Available for Research Groups](https://doc.eresearch.unige.ch/accounting#resources_available_for_research_group) (**HPC account or group member with HPC access required**)

If your research group owns compute nodes but your PI does not appear in the report, please contact the HPC support team.

#### Billing Period

Billing Units are tracked annually from **December 1 to November 30**.

For example:

| Billing Year | Usage Period |
| --- | --- |
| 2026 | 2025-12-01 → 2026-11-30 |
| 2027 | 2026-12-01 → 2027-11-30 |

Monthly reports, annual allocations, remaining balances, and any applicable charges are calculated using this accounting period.
#### Resource Accounting Standardization

All resources are accounted for using a unified Billing Unit (BU) model.
CPU usage, memory consumption, and GPU usage are converted into Billing Units using predefined conversion factors.

For details, please refer to:
[Resource Accounting Uniformization](https://doc.eresearch.unige.ch/accounting#resource_accounting_uniformization)

This unified accounting system allows all resources to be compared and billed consistently, regardless of the underlying hardware.

---
### (Optional) Resource Usage quota

PIs who wish to prevent their research group from exceeding its annual Billing Unit (BU) allocation may request the activation of resource usage limits on their research group account.

When enabled, these limits are enforced by the HPC scheduling system and can help avoid unexpected overconsumption.

Please contact the HPC support team if you would like this option to be enabled for your research group.

---
### Price per Billing Unit

Overview of the SNSF HPC pricing table:

You can download the complete table here, **Please note that the SNSF rate table will be updated soon to reflect the new Billing Unit (BU) terminology**:

University of Geneva users are billed according to category **U1**.

Categories **U2** and **U3** apply to external organizations such as companies and non-academic institutions.

#### Automatic Conversion to Node Pricing
To ensure **fair and predictable pricing**, annual resource consumption is evaluated using two pricing models:

- **Standard BU pricing**
- **Compute node rental pricing**

At the end of the billing period, the **lowest-cost option is automatically applied**.

A standard CPU compute node costs **3,303 CHF/year** and provides approximately **1,342,848 BU/year**.

More details are available in the [Cost of Renting a Compute Node](https://doc.eresearch.unige.ch/hpc/pi_resources#cost_of_renting_a_compute_node) section.

> 
Once annual consumption exceeds **310,382 BU**, billing switches automatically to compute-node pricing.

Compute-node billing is tier-based. Users are billed for the smallest number of compute nodes whose combined annual capacity covers their annual consumption.

| Annual Consumption (BU) | Cost with Standard BU Pricing (U1) | Billing Applied | Annual Cost Charged | Savings | Notes |
| --- | --- | --- | --- | --- | --- |
| **100,000 BU** | **0 CHF** | Included allocation | **0 CHF** | N/A | Covered by the annual UNIGE allocation |
| **150,000 BU** | **785 CHF** | Standard BU pricing | **785 CHF** | **0 CHF** | **100kBU**  Free annual + **50k BU**. Hourly pricing remains more economical |
| **310,382 BU** | **3,303 CHF** | **1 compute node** | **3,303 CHF** | **0 CHF** | Annual consumption threshold triggering automatic compute-node billing |
| **1,342,848 BU** | **21,082 CHF** | **1 compute node** | **3,303 CHF** | **17,779 CHF** | Maximum annual consumption covered by **1 compute node** |
| **1,500,000 BU** | **23,550 CHF** | **2 compute nodes** | **6,606 CHF** | **16,944 CHF** | Consumption exceeds the capacity of **1 compute node** and enters the **2 compute node** tier |
| **2,685,696 BU** | **42,164 CHF** | **2 compute nodes** | **6,606 CHF** | **35,558 CHF** | Maximum annual consumption covered by **2 compute nodes** |
| **4,028,544 BU** | **63,246 CHF** | **3 compute nodes** | **9,909 CHF** | **53,337 CHF** | Maximum annual consumption covered by **3 compute nodes** |
| **5,371,392 BU** | **84,328 CHF** | **4 compute nodes** | **13,212 CHF** | **71,116 CHF** | Maximum annual consumption covered by **4 compute nodes** |
This mechanism can substantially reduce costs for research groups with sustained resource consumption while providing a predictable and transparent billing model.
### Purchasing or Renting Private Compute Nodes

Research groups may purchase or rent compute nodes that are integrated into one of the HPC clusters.

These nodes provide:

- A private partition with higher scheduling priority.
- Reduced queue waiting times.
- Maximum job runtimes of up to **7 days** instead of **4 days** on public resources.
- Additional annual Billing Units associated with the hardware contribution.

Although the hardware is physically attached to a specific cluster, the associated Billing Units can be used across **Baobab**, **Yggdrasil**, and **Bamboo**, including GPU resources.

#### Usage Limits

To preserve the shared nature of the infrastructure, research groups may consume up to **60% of the theoretical annual resource contribution of the hardware they own or rent**.

Example for a compute node with:

- 128 CPU cores
- 512 GB RAM
- 1 year contribution period

**Calculation**:

- **CPU contribution**: 128 × 1.0
- **Memory contribution**: 512 × 0.25
- **Time period**: 24 × 365
- **Usage factor**: 0.6

**Result**:

(128 × 1.0 + 512 × 0.25) × 24 × 365 × 0.6
= **1,342,848 BU**

The resulting allocation (**1,342,848 BU**) can be consumed on any cluster.

The purpose of this model is to ensure fair sharing of the infrastructure while rewarding groups that contribute hardware resources.

#### Rules and Conditions

- **Shared Integration**: The node is integrated into the cluster and remains available to other users when not used by the owning group. See [Partitions](https://doc.eresearch.unige.ch/hpc/slurm#partitions).
- **Usage Policy**: Usage is governed by the 60% allocation rule described above.
- **Administrative Access**: Research groups do not receive administrative access to compute nodes.
- **Maintenance**: Installation and operational maintenance are handled by the HPC team.
- **Decommissioning**: Nodes remain assigned to the associated research group for **five years**. After this period, they are migrated to the public partition and become fully available to the broader HPC community until they reach obsolescence.
- **Warranty** Nodes include a **3-year manufacturer warranty**.
  - After the warranty period, repair costs are the responsibility of the research group.
  - Diagnostic costs charged by the vendor are estimated to **420 CHF**, even if the hardware cannot be repaired.

#### Cost of Renting a Compute Node

The monthly rental cost is calculated from:

- The vendor purchase price.
- A 15% infrastructure surcharge covering racks, cabling, networking, storage, and operational costs.
- A 5-year amortization period.

**Example**:

- Vendor price: 14,361 CHF
- Infrastructure surcharge (15%): 2,154.15 CHF
- Total cost: 16,515.15 CHF
- Amortization period: 60 months

Monthly rental cost:

16,515.15 CHF / 60
= **275.25 CHF** per month

Currently, standard rental nodes consist of:

- 128 AMD CPU cores
- 512 GB RAM

Such a node provides approximately **1.34 million Billing Units per year**.

#### Rental Conditions

- Minimum rental duration: **6 months**
- Unused Billing Units are not carried over and expire at the end of the calendar year.
- Custom hardware configurations may be available upon request.

For quotations or additional information, please contact the HPC support team.

---
### Monthly Usage report for PI

Every PI receives a monthly usage report summarizing the group's resource consumption.

The report includes:

- Annual Billing Units allocated to the research group.
- Billing Units consumed since the beginning of the year.
- Remaining Billing Units.
- Current balance status.

A color indicator helps quickly identify the situation:

- Green: the group remains within its allocated Billing Units.
- Red: the group has exceeded its allocation and additional usage may be billed.

**Example:**

### Monitoring Usage

To monitor your group's current resource consumption, consult:

The dedicated user-friendly web page:

[Baobab OpenOnDemand Dashboard](https://openondemand.baobab.hpc.unige.ch/pun/sys/ug_my_hpc_usage/) (**HPC account or group member with HPC access required**)

or  with linux CLI:

[Report and Statistics with sreport](https://doc.eresearch.unige.ch/accounting#report_and_statistics_with_sreport) (**HPC account or group member with HPC access required**)
