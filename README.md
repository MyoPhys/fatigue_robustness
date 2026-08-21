# fatigue_robustness
August 20, 2026

This repository contains the raw data files and R scripts necessary to replicate the results of the study
*The storage and release of energy from elastic tissues is unaffected by muscle fatigue* by Jeffrey P. Olberding. This is a study of the effects of muscle fatigue on the ability of elastic tissues to store and release energy during isotonic and Latch-mediated Spring Actuation (LAMSA) type contractions.

**Data use license agreement: These data are provided only for reviewers and readers to check and replicate my analyses. For use of the data in other original research, please contact Jeffrey Olberding for permission.**

Analyses were performed using RStudio V.2024.04 (Posit Software, PBC) running R V.4.4.1 (R Core Team).

**1. "Raw Data"** - folder containing raw data files from contractile experiments.
   Data files are names using the convention Date_Muscle number_Type of contraction_Before or after fatigue.
   
**2. "Metadata.csv"** - csv file containing data unique to each recording needed for analysis. The following columns are defined:
   - Date - the date the data were collected
   - Indiv - The individual identifiers of the frogs used for data collection
   - Muscle - Indicating the 1st and 2nd muscle used from each individual
   - Load - Code indicating the type of contraction and the load used. E.g. "Iso20" is an isotonic contraction set to 20% of maximum isometric tension.
   - Before/After - Coding whether the recorded data were before or after the fatigue protocol.
   - Muscle Mass (kg) - The mass of the muscle used for recording
   - MTU Length (mm) - The length of the muscle used for recording
   - Force Threshold (N) - The % of maximum isometric tension used to define the load of the contraction.
**3. "Analyze ISO Fatigue V3.r"** - Script used to analyze an individual isotonic type contraction.
   - Raw data files of isotonic contractions are indicated by a file name containing "Iso", e.g. "*10_02_23_Muscle1_Iso20_After.txt*"
   - This script requires edits to include metadata matching the individual trials being analyzed. Specifically, MTU Length from Metadata.csv must be entered in line 6.
   - Outputs the following variables: averageforce, tendonlength, maxshortening, maxcontraction, timetoshort, loadenergy, peakfasciclev, peakfasciclep
   - All raw files for isotonic type contractions were analyzed and data were compiled into *Isotonic Data V2.csv*.
   - This script will create a plot of the data for the individual trial being run and was used in the creation of Figure 3.
     
**4. "Analyze LAMSA Fatigue V3.r"** - Script used to analyze an individual LaMSA type contraction.
   - Raw data tiles of LAMSA type contractions are indicated by a file name containing "LAMSA", e.g. "*10_02_23_Muscle1_LAMSA20_After.txt*"
   - This script requires edits to include metadata matching the individual trials being analyzed. Specifically, MTU Length, Muscle Mass, and Force Threshold from Metadata.csv must be entered in lines 6, 7, and 8, respectively.
   - Outputs the following variables: max_tendon, max_force, load_pre, overshoot, e_return,direct,  peak_tendonp, peak_fasciclep, peak_MTUp, peak_fasciclev, loading_power, loading_time, recoiltime
   - All raw files for LAMSA type contractions were analyzed and data were compiled into *LAMSA Data V3.csv*
   - This script will create a plot of the data for the individual trial being run and was used in the creation of Figures 4 and 5.
     
**5. "Isotonic Data V2.csv"** - File containing the compiled output from analyzing raw data files using *Analyze ISO Fatigue V3.r*.
     - Variables include metadata, as well as, averageforce, tendonlength, maxshortening, maxcontraction, timetoshort, loadenergy, peakfasciclev, peakfasciclep output from *Analyze ISO Fatigue V3.r* and tendon stiffness, mms elastic, mms muscle, mms p calculated from those output variables.
     
**6. "Analyze Isotonic Results v2.r"** - Script use to generate descriptive statistics, run ANOVAs, and create box plots the compiled isotonic data.
  - This script runs file *Isotonic Data V2.csv*
   - This script outputs descriptive statistics for variables. Lines 24-31 must be edited by replacing "tendonp" with the variable of interest. These were used to create Table 2.
     - This script runs ANOVAs on all variables and outputs statistical summaries.
     - This script creates boxplots for all variables, which were used to create Figure 6.

**7. "LAMSA Data V3.csv"** - File containing the compiled output from analyzing raw data files using *Analyze LAMSA Fatigue V3.r*.
     - Variables include metadata, as well as, max_tendon, max_force, load_pre, overshoot, e_return,direct,  peak_tendonp, peak_fasciclep, peak_MTUp, peak_fasciclev, loading_power, loading_time, recoiltime output from *Analyze LAMSA Fatigue V3.r* and mms.pre, mms.over, mms.return, mms.direct, tendon.stiffness, efficiency calculated from those output variables.

**8. "Analyze LAMSA Results v2.r"** - Script use to generate descriptive statistics, run ANOVAs, and create box plots the compiled LAMSA data.
  - This script runs file *LAMSA Data V3.csv*
   - This script outputs descriptive statistics for variables. Lines 24-31 must be edited by replacing "tendonp" with the variable of interest. These were used to create Table 4.
     - This script runs ANOVAs on all variables and outputs statistical summaries.
     - This script creates boxplots for all variables, which were used to create Figure 6.
