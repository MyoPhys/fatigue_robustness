# fatigue_robustness
Files and code for original manuscript on fatigue robustness

1. "Raw Data" - folder containing raw data files from contractile experiments.
   Data files are names using the convention Date_Muscle number_Type of contraction_Before or after fatigue.
   
2. "Metadata.csv" - csv file containing data unique to each recording needed for analysis. The following columns are defined:
   - Date - the date the data were collected
   - Indiv - The individual identifiers of the frogs used for data collection
   - Muscle - Indicating the 1st and 2nd muscle used from each individual
   - Load - Code indicating the type of contraction and the load used. E.g. "Iso20" is an isotonic contraction set to 20% of maximum isometric tension.
   - Before/After - Coding whether the recorded data were before or after the fatigue protocol.
   - Muscle Mass (kg) - The mass of the muscle used for recording
   - MTU Length (mm) - The length of the muscle used for recording
   - Force Threshold (N) - The % of maximum isometric tension used to define the load of the contraction.
3. "Analyze ISO Fatigue V3.r" - Script used to analyze an individual isotonic type contraction.
   - Raw data files of isotonic contractions are indicated by a file name containing "Iso", e.g. "*10_02_23_Muscle1_Iso20_After.txt*"
   - This script requires edits to include metadata matching the individual trials being analyzed. Specifically, MTU Length from Metadata.csv must be entered in line 6.
   - Outputs the following variables: averageforce, tendonlength, maxshortening, maxcontraction, timetoshort, loadenergy, peakfasciclev, peakfasciclep
   - All raw files for isotonic type contractions were analyzed and data were compiled into *Isotonic Data V3.csv*.
   - This script will create a plot of the data for the individual trial being run and was used in the creation of Figure 3.
6. "Analyze LAMSA Fatigue V3.r" - Script used to analyze an individual LaMSA type contraction.
   - Raw data tiles of LAMSA type contractions are indicated by a file name containing "LAMSA", e.g. "*10_02_23_Muscle1_LAMSA20_After.txt*"
   - This script requires edits to include metadata matching the individual trials being analyzed. Specifically, MTU Length, Muscle Mass, and Force Threshold from Metadata.csv must be entered in lines 6, 7, and 8, respectively.
   - Outputs the following variables: max_tendon, max_force, load_pre, overshoot, e_return,direct,  peak_tendonp, peak_fasciclep, peak_MTUp, peak_fasciclev, loading_power, loading_time, recoiltime
   - All raw files for LAMSA type contractions were analyzed and data were compiled into *LAMSA Data V3.csv*
   - Data were compiled in creation of Table 3. Individual plots were used in the creation of Figures 4 and 5.
8. "Analyze ISO Results v2.r" - code to run ANOVAs and create boxplots for isotonic data. Used to create Table 2 and Figure 6.
9. "Analyze LAMSA Results v2.r" - code to run ANOVAs and create boxplots for LaMSA data. Used to create Table 4 Figure 6.
