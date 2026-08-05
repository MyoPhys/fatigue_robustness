# fatigue_robustness
Files and code for original manuscript on fatigue robustness

1. "Raw Data" - folder containing raw data files from contractile experiments
2. "Metadata.csv" - csv file containing data unique to each recording needed for analysis. The following columns are defined:
   a. Date - the date the data were collected
   b. Indiv - The individual identifiers of the frogs used for data collection
   c. Muscle - Indicating the 1st and 2nd muscle used from each individual
   d. Load - Code indicating the type of contraction and the load used. E.g. "Iso20" is an isotonic contraction set to 20% of maximum isometric tension.
   e. Before/After - Coding whether the recorded data were before or after the fatigue protocol.
   f. Muscle Mass (kg) - The mass of the muscle used for recording
   g. MTU Length (mm) - The length of the muscle used for recording
   h. Force Threshold (N) - The % of maximum isometric tension used to define the load of the contraction.
4. "Analyze ISO Fatigue V3.r" - code to analyze an individual isotonic type contraction; also creates plots of single contraction data. Requires edits for metadata matching the individual trials being analyzed. MTU Length from Metadata.csv must be entered in line 6. Data were compiled in creation of Table 1. Individual plots were used in the creation of Figure 3.
5. "Analyze LAMSA Fatigue V3.r" - code to analyze and individual LaMSA type contraction; also creates plots of single contraction data. Requires edits for metadata matching the individual trials being analyzed. MTU Length, Muscle Mass, and Force Threshold from Metadata.csv must be entered in lines 6, 7, and 8, respectively. Data were compiled in creation of Table 3. Individual plots were used in the creation of Figures 4 and 5.
6. "Analyze ISO Results v2.r" - code to run ANOVAs and create boxplots for isotonic data. Used to create Table 2 and Figure 6.
7. "Analyze LAMSA Results v2.r" - code to run ANOVAs and create boxplots for LaMSA data. Used to create Table 4 Figure 6.
