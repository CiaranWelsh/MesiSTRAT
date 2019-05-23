Modelling meeting
=================
- Problem with FCS
    - 5% to 2%
    - Cells kept at low FCS, stimulate for short period of time to see the dynamics.
    - Problem is that in full serum, we loose the dynamics.
    - INF-gamma
    - How long do you need to do these experiments for?
        - Measure metabolites, time needed for acucmulationi.
        - Singalling = faster, metabolites = slower
        - Idea: increase cell concentration to increase metabolites
        - metabolites
            - turnover and concentration are different. YOu may not see changesin concentration but you may have differences in turnover.
- Short term, we need to get some first models parameterised - end of june.
- Stary working with method that we have an works,
- How do we model growth factors?
    - Assume constants or measure
- Tryprophan, metabolites of tryprophan, growth factors?
- Glucose
    - too high conc. decrease would be small and we do not get reliable flux
    - Depends on time scale. Long - this becomes important. Short, not so important
- if you have ido1 induction, tryptophan will run out eventualy.
- Tryptophan deprevation is important in these tumours.
    - How much deprevation is physiological ?
    - Decision, measure trypt. keep this constant.
        - media, cell count, culture dishes all will affect the tryptophan output and must be kept constant
        - This is a reference condition. Modulate conditions and measure how this impacts the output of tryptophan
        - MEasure these metabolites in the medium.


- So tryptophan is important in BC tumours. This is what we are interested in evaluating.
    - How do we measure the importance of tryptophan?
    - What can we use as a proxy for the effect of tryptophan on tumours?

- 25mmol glucose
- Stay at 2%FCS
- Time courses run over experiments.
    - measure status of singlaling and metabolites
    - we expecte that once were belfow the amount of trypt we might also see shifts in other parts othe network
- IGF EGF, combination
- Fulvestrant
- mtor/MAPK pathwyas
- AMPK
- glucose, trypt+metabolites,
- Output
- MCF7 we see an increase in TOR signalling when trypt is low
    - And Mek Erk.
- Se were using mTOR/Erk activation as a proxy for tumour progression.
- metabolites experiments
- Literature data, alternative data and conditions that are useful for informing the model.
- To proteomics or not to proteomics.
- Month 24, ralative quantification TMT based quantification to look at proteins in our networks
- Tarteded TMT proteomics
- identifyicatoin of subgroups
    - Setup models basedon different cell lines
    - They behave differently
    - SA. Which are most sensitive
    - Define set of minimal [arameters to be capable of transferring data/models between patients
    - How do you derive the patients
        - based on cell lines?
        - You need to do some kind of clustering
        - signatures of different cells
    - Initial = bottom up. Later we do clustering. But when? Afterwards
- First test case. Patient derived cell lines. Different response to endocrine therapy (cultures xenografts).
- Motility could be a readout
- When we run the cells under these condition we do not see much
- They do see an effect on the motility of these cells
- Drugs that target ER have an effect on motility.
- What defines the subgroups?
    - Find different dynamics in th emodels
    - Confirm n PDx's
    - How many patients are required to sub group them?
    - MAPK/mTOR are directly targets by the drugs we are using.
    - TRy/NAD are improtant for drug response because they are modulatng these oncogenic signalling netowkr
    - MEasuring components of these networks may help us stratify these patients.
    - Factors whos expression will not alter network behaviour.
    - Complicated of feature engineering.
- Linking signalling with metabolic. Interventions are needed.
- We want to know which markers enable us to discriminate between patient subtypes.





important points
----------------
- The key issues are:
    - What do we stimulate with
        - conditions, which, how many?
        - Fulvestrant, INF-y
        -
    - What time points do we measure at?
    - which are the most important nodes in the network
    - What is the model output?
    - Network is known.
        - Get hold of this and build the network
    -
- Output?
- What is the point?
- Tumour size? IDO1? Tryptophan?




Summary of plan
---------------
- First model by middle of 2019. 
- 





























