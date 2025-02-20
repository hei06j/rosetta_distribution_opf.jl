# NB readme copied from full dataset, all data available from
https://doi.org/10.25919/jaae-vc35


# Three-wire and Four-wire variants of ENWL feeders

This dataset contains three-wire and four-wire variants of the ENWL low-voltage networks.


## Original dataset

The original dataset consists of 128 low-voltage networks, complete with load profiles. For a detailed discussion, refer to ENWL - Low Voltage Network Solutions project [website](https://www.enwl.co.uk/future-energy/innovation/smaller-projects/low-carbon-networks-fund/low-voltage-network-solutions/).


For the original networks, the line impedance is specified in sequence components. This implies that two simplifications have been applied:
1. a Kron-reduction, which assumes that the neutral wire is grounded everywhere;
2. discarding off-diagonal elements of the series impedance matrix in sequence coordinates.

This dataset however, contains three-wire and four-wire variants of the original feeders. The simplifications which we discussed above, are irreversible. Therefore, we developed a procedure to re-instantiate the dataset with reasonable detailed line models. 



## Four-wire extension

First, we compiled a library of detailed line models. These line models were provided by [1]. Next, for each line in the datatset, we assign the four-wire line model which matches the original one most closely in terms of positive sequence impedance [2]. Further numerical experiments were performed in [3].


## Three-wire variants

We next transform the four-wire line models to the following three-wire models and make necessary changes to the feeder models.

### Kron-reduced line models

For more details on this line model, please refer to [2,3].

### Phase-neutral line models

For more details on this line model, please refer to [3].

### Modified phase-neutral line models

For more details on this line model, please refer to [3].


## Usage

This dataset contains the following folders:
- `four-wire`: These are the four-wire extensions of each feeder in the dataset with four-wire line models.
- `three-wire-Kron-reduced`: These are the three-wire extensions of each feeder in the dataset with kron-reduced line models.
- `three-wire-phase-to-neutral`: These are the three-wire extensions of each feeder in the dataset with phase-neutral line models.
- `three-wire-modified-phase-to-neutral`: These are the three-wire extensions of each feeder in the dataset with modified phase-neutral line models.


Note that the loads default to 1 kW in the output file, and the dss loadshape files are not included. You can refer to the original dataset [website](https://www.enwl.co.uk/future-energy/innovation/smaller-projects/low-carbon-networks-fund/low-voltage-network-solutions/) to obtain load profiles for various scenarios.


## References

[1] Urquhart, Andrew J., and Murray Thomson. 2019. “Cable Impedance Data”. figshare. https://hdl.handle.net/2134/15544. Creative Commons Attribution-NonCommercial 3.0 Unported License.

[2] 'Distribution Network Modeling: From Simulation Towards Optimization, Sander Claeys, Phd Dissertation, KU Leuven, Belgium 2021.

[3] Frederik Geth, Rahmat Heidari, and Arpan Koirala. 2022. Computational analysis of impedance transformations for four-wire power networks with sparse neutral grounding. In Proceedings of the Thirteenth ACM International Conference on Future Energy Systems (e-Energy '22). Association for Computing Machinery, New York, NY, USA, 105–113. https://doi.org/10.1145/3538637.3538844