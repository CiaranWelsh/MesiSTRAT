import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import os, glob, pandas, numpy
try:
    from pycotools import model, viz, tasks
except ImportError:
    print('Warning: Pycotools is not installed. Are you running Python3? ')

import tellurium as te
import seaborn
from collections import OrderedDict
from itertools import cycle
from decimal import Decimal, getcontext
from multiprocessing import Process, Pool, cpu_count
from collections import Counter


"""
Idea: Simulate 100 models with random parameters. Record whether each model satisifies each condition. Or better 
yet, compute a semi RSS. Then use the values of the parameters to as a regressor neural network for prediction of 
best parameters for satisfying certain conditions. 
"""

seaborn.set_style('white')
seaborn.set_context('talk', font_scale=1)


MODEL_CODE                      = 'E_A_48'
CURRENT_SPECIES = [
    'pSmad2Tot', 'pSmad2', 'pSmad2n', 'Smad2n', 'Smad2Tot',
    'pErk', 'ppErk', 'pAkt', 'pS6K', 'pmTORC1',
]
CURRENT_SPECIES = [
    'pSmad2Tot', 'pMek', 'ppMek',
    'pErk', 'ppErk', 'pAkt', 'pS6K', 'pmTORC1',
]
# CURRENT_SPECIES = ['pSmad2Tot']
# CURRENT_SPECIES = [
#     'pSmad2Tot', 'pSmad2', 'pSmad2n', 'Smad2n', 'Smad2Tot',
#
# ]

SIMULATE_TIME_SERIES            = False
SIMULATE_BAR_GRAPHS             = False
OPEN_CONDITION_WITH_COPASI      = False

QUALITATIVE_FITTING             = True
OPTIMIZE                        = False

GET_PARAMETERS_FROM_COPASI      = False

CONFIGURE_PARAMETER_ESTIMATION  = False

DOSE_RESPONSE_GROWTH_FACTOR     = False
DOSE_RESPONSE_TGFB              = False
GET_ODES_WITH_ANTIMONY          = False
GET_MODEL_AS_SBML               = False
SIMULATE_INPUTS                 = False

"""
These are arguments for the conditions simulation functions. 
    :param GF: starting amount of GrowthFactors
    :param pretreatment: either 'AZD', or 'MK2206'. This is added as an event at the time specified by pretreatment time
    :param pretreatment_time: The time at which the variable specified by 'pretreatment' is added
    :param EV: Starting amount of Everolimus
    :param serum_starve_event: Boolean, whether to remove serum, aka GrowthFactors by event
    :param TGFb_event: Boolean. Whether to add 1 unit TGF at t=-45min (71.75h)

"""

AZD_CONDITIONS = OrderedDict()
AZD_CONDITIONS['D']         =   [1,  None,     None,       0,  True,   False]
AZD_CONDITIONS['T']         =   [1,  None,     None,       0,  True,   True]
AZD_CONDITIONS['E']         =   [1,  None,     None,       1,  True,   True]
AZD_CONDITIONS['E_A_72']    =   [1,  "AZD",    0.0,        1,  True,   True]
AZD_CONDITIONS['E_A_48']    =   [1,  "AZD",    24,         1,  True,   True]
AZD_CONDITIONS['E_A_24']    =   [1,  "AZD",    48.0,       1,  True,   True]
AZD_CONDITIONS['E_A_1.25']  =   [1,  "AZD",    70.75,      1,  True,   True]
AZD_CONDITIONS['A_72']      =   [1,  "AZD",    0.0,        0,  True,   True]
AZD_CONDITIONS['A_48']      =   [1,  "AZD",    24,         0,  True,   True]
AZD_CONDITIONS['A_24']      =   [1,  "AZD",    48.0,       0,  True,   True]
AZD_CONDITIONS['A_1.25']    =   [1,  "AZD",    70.75,      0,  True,   True]

MK_CONDITIONS = OrderedDict()
MK_CONDITIONS['D']          =     [1, None,      None,   0, True, False]
MK_CONDITIONS['T']          =     [1, None,      None,   0, True, True]
MK_CONDITIONS['E']          =     [1, None,      None,   1, True, True]
MK_CONDITIONS['E_M_72']     =     [1, "MK2206",  0.0,    1, True, True]
MK_CONDITIONS['E_M_48']     =     [1, "MK2206",  24.0,   1, True, True]
MK_CONDITIONS['E_M_24']     =     [1, "MK2206",  48.0,   1, True, True]
MK_CONDITIONS['E_M_1.25']   =     [1, "MK2206",  70.75,  1, True, True]
MK_CONDITIONS['M_72']       =     [1, "MK2206",  0.0,    0, True, True]
MK_CONDITIONS['M_48']       =     [1, "MK2206",  24.0,   0, True, True]
MK_CONDITIONS['M_24']       =     [1, "MK2206",  48.0,   0, True, True]
MK_CONDITIONS['M_1.25']     =     [1, "MK2206",  70.75,  0, True, True]


MODEL_SPECIES = ['TGFbR',   'TGFbR_a',  'TGFbR_EE', 'TGFbR_Cav',
                 'Smad2',   'pSmad2',   'Raf',      'pRaf',
                 'Mek',     'pMek',     'ppMek',    'Erk',        'pErk',
                 'ppErk',   'PI3K',     'pPI3K',    'Akt',        'pAkt',
                 'mTORC1',  'pmTORC1',  'S6K',      'pS6K',
                 'pSmad2Tot', 'pSmad2n', 'Smad2Tot', 'Smad2n'
                 ]

MODEL_INPUTS = ['TGFb', 'Everolimus', 'MK2206', 'AZD', 'GrowthFactors']

MODEL_SPECIES += MODEL_INPUTS

ALL_CONDITIONS = ['D', 'T', 'E', 'E_A_72', 'E_A_48',
                  'E_A_24', 'E_A_1.25', 'A_72', 'A_48',
                  'A_24', 'A_1.25', 'E_M_72',
                  'E_M_48', 'E_M_24', 'E_M_1.25', 'M_72',
                  'M_48', 'M_24', 'M_1.25']



def cross_talk_model_antstr():
    """
    Antimony string model containing the most currect version of
    the CrossTalk model from Katrine and Patricia's project.
    :return:
    """
    return """
    function MM(km, Vmax, S)
        Vmax * S / (km + S)
    end

    function MMWithKcat(km, kcat, S, E)
        kcat * E * S / (km + S)
    end
    

    function NonCompetitiveInhibition(km, ki, Vmax, n, I, S)
        Vmax * S / ( (km + S) * (1 + (I / ki)^n ) )
    end
    
    function MA1(k, S)
        k * S
    end
    
    function MA2(k, S1, S2)
        k * S1 * S2
    end
    
    function MA1Mod(k, S, M)
        k * S * M
    end
    
    function MA2Mod(k, S1, S2, M)
        k * S1 * S2 * M
    end
    
    function CompetitiveInhibition(km, ki, kcat, E, I, S)
        kcat * E * S / (km + S + ((km * I )/ ki)  )
    end
    
    model TGFbModule()
        compartment Cell = 1.0
        
        var TGFbR           in Cell  
        var TGFbR_a         in Cell  
        var TGFbR_EE        in Cell
        var TGFbR_Cav       in Cell
        var Smad2           in Cell  
        var pSmad2          in Cell  
        var pSmad2n         in Cell
        var Smad2n          in Cell 
        var Mek             in Cell
        var pMek            in Cell  
        var Erk             in Cell
        var pErk            in Cell  
        var PI3K            in Cell  
        var pPI3K           in Cell  
        var Akt             in Cell
        var pAkt            in Cell  
        var mTORC1          in Cell  
        var pmTORC1         in Cell  
        var S6K             in Cell
        var pS6K            in Cell  
        var Raf             in Cell
        var pRaf            in Cell
        var pMek            in Cell
        var ppMek           in Cell
        var pMek            in Cell
        var pErk            in Cell
        var ppErk           in Cell
        
        
        const TGFb             in Cell
        const AZD              in Cell
        const GrowthFactors    in Cell
        const MK2206           in Cell
        const Everolimus       in Cell
    
    
        //TGFb module
        TGF_R1 : TGFbR       => TGFbR_a    ; Cell * kTGFbOn       *TGFbR      *TGFb      ;
        TGF_R2 : TGFbR_a     => TGFbR      ; Cell * kTGFbOff      *TGFbR_a               ;
        TGF_R3 : TGFbR_a     => TGFbR_EE   ; Cell * kTGFbRIntern  *TGFbR_a               ; 
        TGF_R4 : TGFbR_EE    => TGFbR_a    ; Cell * kTGFbRRecyc   *TGFbR_EE              ;
        TGF_R5 : TGFbR_a     => TGFbR_Cav  ; Cell * kTGFbRIntern  *TGFbR_a               ;
        TGF_R6 : TGFbR_Cav   => TGFbR_a    ; Cell * kTGFbRRecyc   *TGFbR_Cav             ;
        TGF_R8 : Smad2       => pSmad2     ; Cell * MMWithKcat(kSmad2Phos_km, kSmad2Phos_kcat, Smad2, TGFbR_EE );
        TGF_R10: pSmad2      => pSmad2n    ; Cell  * MM(kpSmad2Imp_km, kpSmad2Imp_Vmax, pSmad2                  );
        TGF_R11: pSmad2n     => pSmad2     ; Cell * MM(kSmad2Exp_km, kSmad2Exp_Vmax, pSmad2n                 );
        TGF_R12: Smad2       => Smad2n     ; Cell  * MM(kSmad2Imp_km, kSmad2Imp_Vmax, Smad2                     );
        TGF_R13: Smad2n      => Smad2      ; Cell * MM(kSmad2Exp_km, kSmad2Exp_Vmax, Smad2n                    );
        TGF_R14: pSmad2n     => Smad2n     ; Cell  * MM(kpSmad2Dephos_km, kpSmad2Dephos_Vmax, pSmad2n);
        
        //MAPK module
        MAPK_R0  : Raf     => pRaf      ; Cell*GrowthFactors*NonCompetitiveInhibition(kRafPhos_km,  kRafPhos_ki, kRafPhos_Vmax, kRafPhos_n, ppErk, Raf);
        MAPK_R1  : pRaf    => Raf       ; Cell*MM(            kRafDephos_km ,   kRafDephosVmax,      pRaf           );
        MAPK_R2  : Mek     => pMek      ; Cell*CompetitiveInhibition(    kMekPhos_km1 , kMekPhos_ki1, kMekPhos_kcat1, pRaf, AZD, Mek       );
        MAPK_R3  : pMek    => ppMek     ; Cell*CompetitiveInhibition(    kMekPhos_km1 , kMekPhos_ki1, kMekPhos_kcat1, pRaf, AZD, pMek     );
        MAPK_R4  : ppMek   => pMek      ; Cell*MM(            kMekDephos_km1,   kMekDephos_Vmax1,     ppMek         );
        MAPK_R5  : pMek    => Mek       ; Cell*MM(            kMekDephos_km1,   kMekDephos_Vmax1,     pMek          );
        MAPK_R6  : Erk     => pErk      ; Cell*MMWithKcat(    kErkPhos_km1,     kErkPhos_kcat1, Erk,  ppMek         );
        MAPK_R7  : pErk    => ppErk     ; Cell*MMWithKcat(    kErkPhos_km1,     kErkPhos_kcat1, pErk, ppMek         );
        MAPK_R8  : ppErk   => pErk      ; Cell*MM(            kErkDephos_km1,   kErkDephos_Vmax1,     ppErk         );
        MAPK_R9  : pErk    => Erk       ; Cell*MM(            kErkDephos_km1,   kErkDephos_Vmax1,     pErk          );

        
        //PI3K Module
        //n(km, ki, kcat, E, I, S)
        PI3K_R1 :   PI3K    => pPI3K    ;  Cell *  kPI3KPhosByGF       *PI3K       *GrowthFactors ;
        PI3K_R2 :   pPI3K   => PI3K     ;   Cell *  kPI3KDephosByS6K    *pPI3K      *pS6K        ;
        PI3K_R3 :   Akt    => pAkt      ;   Cell *  CompetitiveInhibition(kAktPhos_km, kAktPhos_ki, kAktPhos_kcat, pPI3K, MK2206, Akt)              ;
        PI3K_R4 :   pAkt    => Akt      ;   Cell *  MM(kAktDephos_km, kAktDephos_Vmax,pAkt)         ;
        PI3K_R5_1 :   mTORC1 => pmTORC1   ;   Cell *  CompetitiveInhibition(kmTORC1Phos_km, kmTORC1Phos_ki, kmTORC1Phos_kcat, pAkt, Everolimus, mTORC1)  ;
        //PI3K_R5_2 :   mTORC1 => pmTORC1   ;   Cell *  kmTORC1PhosByTSCs * mTORC1  ;
        PI3K_R6 :   pmTORC1 => mTORC1   ;   Cell *  MM(kmTORC1Dephos_km, kmTORC1Dephos_Vmax, pmTORC1);
        PI3K_R7 :   S6K     => pS6K     ;   Cell *  MMWithKcat(kS6KPhosBymTORC1_km, kS6KPhosBymTORC1_kcat, S6K, pmTORC1) ;
        PI3K_R8 :   pS6K    => S6K      ;   Cell *  MM(kS6KDephos_km,kS6KDephos_Vmax, pS6K)                    ;
        
        // Cross talk reactions
        CrossTalkR1  :    Raf       => pRaf      ;   Cell *  MMWithKcat(kRafPhosByTGFbR_km, kRafPhosByTGFbR_kcat, Raf, TGFbR_Cav)    ;
        CrossTalkR2  :    Raf       => pRaf      ;   Cell *  MMWithKcat(kRafPhosByPI3K_km,kRafPhosByPI3K_kcat, Raf, pPI3K)           ;
        CrossTalkR3  :    PI3K      => pPI3K     ;   Cell *  MMWithKcat(kPI3KPhosByTGFbR_km, kPI3KPhosByTGFbR_kcat, PI3K, TGFbR_Cav) ;
        CrossTalkR4  :    pPI3K     => PI3K      ;   Cell *  kPI3KDephosByErk    *pPI3K      *ppErk        ;
        CrossTalkR5  :    Smad2     => pSmad2    ;   Cell *  MMWithKcat(kSmad2PhosByAkt_km, kSmad2PhosByAkt_kcat, Smad2, pAkt)       ;
        CrossTalkR6  :    pSmad2n   =>  Smad2n   ;   Cell *  MMWithKcat(kSmad2DephosByErk_km, kSmad2DephosByErk_kcat, pSmad2n, ppErk)       ;
        
        // Species initializations:
        TGFbR       = 76.8396790634687;
        TGFbR_a     = 0.966718661034664;
        TGFbR_EE    = 12.55032566215;
        TGFbR_Cav   = 9.66718661034664;
        Smad2       = 49.6779
        pSmad2      = 0.500563;
        pSmad2n     = 0.0447836;
        Smad2n      = 49.7768;
        Smad2Tot    := pSmad2 + Smad2 + pSmad2n + Smad2n
        pSmad2Tot   := pSmad2 + pSmad2n;
        Mek         = 252.876273823102;
        pMek        = 56.0642557607438;
        Erk         = 183.331905604514;
        pErk        = 117.882228538386;
        PI3K        = 99.9603027531632;
        pPI3K       = 0.039697246836694;
        Akt         = 99.413240815732;
        pAkt        = 0.586759181268058;
        mTORC1      = 90.2724562462086;
        pmTORC1     = 9.72754375379144;
        S6K         = 90.5069860202662;
        pS6K        = 9.49301398073377;
        Raf         = 87.932595095561;
        pRaf        = 12.0674049063339;
        ppMek       = 29.0295122061532;
        ppErk       = 43.0862638108839;
        
        // Variable initializations:                    
        TGFb                    = 0.005;                    
        GrowthFactors           = 1;                    
        kTGFbOn                 = 0.1;                  
        kTGFbOff                = 0.04;                 
        kTGFbRIntern            = 0.3333333333;                 
        kTGFbRRecyc             = 0.03333333333;                    
        
        kSmad2Phos_km           = 50;                  
        kSmad2Phos_kcat         = 0.1;                    
        kSmad2PhosByAkt_km      = 40;                  
        kSmad2PhosByAkt_kcat    = 0.1;                   
        kpSmad2Dephos_km         = 60;                 
        kpSmad2Dephos_Vmax       = 65;
        kSmad2DephosByErk_km    = 30;                   
        kSmad2DephosByErk_kcat  = 7.5;   
        
        mul = 3
        kSmad2Imp_km           = 90;
        kpSmad2Imp_km         := kSmad2Imp_km 
        kpSmad2Imp_Vmax       := kSmad2Imp_Vmax * mul
        kSmad2Imp_Vmax         = 38.466;
        kSmad2Exp_km           = 20;
        kSmad2Exp_Vmax         = 20;
        
        kRafPhos_km             = 10;                   
        kRafPhos_ki             = 3.5;                  
        kRafPhos_Vmax           = 9000;                 
        kRafPhos_n              = 1;                    
        kRafDephos_km           = 8;                    
        kRafDephosVmax          = 3602.5;                   
        kMekPhos_km1            = 15;                   
        kMekPhos_ki1            = 0.25;  //original 0.25                 
        kMekPhos_kcat1          = 90;                  
        AZD                     = 0;                    
        kMekDephos_km1          = 15;                   
        kMekDephos_Vmax1        = 2700;                 
        kErkPhos_km1            = 50;                  
        kErkPhos_kcat1          = 200;                 
        kErkDephos_km1          = 15;                   
        kErkDephos_Vmax1        = 1800;                 
        kPI3KPhosByGF           = 0.239474698704283;                    
        kPI3KDephosByS6K        = 25;                   
        kAktPhos_km             = 12.5;                 
        kAktPhos_ki             = 0.01;                 
        kAktPhos_kcat           = 1.5; // original: 2.9215;                 
        MK2206                  = 0;                    
        kAktDephos_km           = 15;                   
        kAktDephos_Vmax         = 30;                   
        kmTORC1Phos_km          = 3;                    
        kmTORC1Phos_ki          = 0.001;                    
        kmTORC1Phos_kcat        = 0.35;        
        Everolimus              = 0;                    
        kmTORC1Dephos_km        = 100;                  
        kmTORC1Dephos_Vmax      = 1;                    
        kS6KPhosBymTORC1_km     = 100;                  
        kS6KPhosBymTORC1_kcat   = 0.5;                  
        kS6KDephos_km           = 10;                   
        kS6KDephos_Vmax         = 50;                   
        kRafPhosByTGFbR_km      = 25;                   
        kRafPhosByTGFbR_kcat    = 265;                  
        kRafPhosByPI3K_km       = 50;                   
        kRafPhosByPI3K_kcat     = 5;                    
        kPI3KPhosByTGFbR_km     = 10;                   
        kPI3KPhosByTGFbR_kcat   = 50;                   
        kPI3KDephosByErk        = 0.5;                  
        //kTGFbRInternByAkt     = 0.01;           

        unit volume = 1 litre;
        unit time_unit = 3600 second;
        unit substance = 1e-9 mole;
        
end
"""


def load_model_with_pyco(ant, copasi_filename):
    """
    Load an antimony model into copasi using Pycotools wrapper around
    tellurium functions
    :param copasi_filename:
    :param ant:
    :return:
    """
    with model.BuildAntimony(copasi_filename) as loader:
        mod = loader.load(
            ant
        )
    return mod


def dose_response(model_string,
                  scan_parameter,
                  start_value,
                  end_value,
                  number_of_points,
                  selection,
                  log10=False
                  ):
    """
    Compute a dose response curve for the antimony model in `model_string`.
    :param model_string: Antimony model
    :param scan_parameter: Which model component to use as scan variable
    :param start_value: Starting value of scan
    :param end_value: End value of scan
    :param number_of_points: Number of intervals
    :param selection: Which model components to track. This is a list of model components. Strings are automatically
    converted into a 1 element list.
    :return:
    """
    mod = te.loada(model_string)
    mod.conservedMoietyAnalysis = True
    if not isinstance(selection, list):
        selection = [selection]

    p = te.SteadyStateScan(mod,
                           value=scan_parameter,
                           startValue=start_value,
                           endValue=end_value,
                           numberOfPoints=number_of_points,
                           selection=selection
                           )
    res = p.steadyStateSim()
    if log10:
        res[:, 0] = numpy.log10(res[:, 0])

    # print(res)
    fig = plt.figure()
    for i in range(len(selection)):
        plt.plot(res[:, 0], res[:, i + 1], label=selection[i], marker='o')
    seaborn.despine(fig=fig, top=True, right=True)
    plt.xlabel('Initial Conc. of {}'.format(scan_parameter))
    plt.ylabel(r'$nmol L^{-1}$')
    plt.legend(loc='best')
    return fig


def __dep__dose_response_for_growth_factors():
    """
    Now deprecated since GrowthFactors is now boolean
    Use the dose repsonse function to see the effect of changing the
    initial concentration of GrowthFactors. Plot some output
    :return:
    """
    l = [
        ['PI3K', 'pPI3K'],
        ['Mek', 'pMek'],
        ['Erk', 'pErk'],
        ['Akt', 'pAkt'],
        ['S6K', 'pS6K'],
    ]

    for i, j in l:
        # print(mod.model.keys())
        fig = dose_response(cross_talk_model_antstr(), 'GrowthFactors', 1, 10000, 1000, [i, j])
        plt.title('Dose Response of {} and \n{} to GrowthFactors (log10 nM)'.format(i, j))
        fname = os.path.join(graphs_directory, 'GrowthFactorsDoseResponse{}.png'.format(i))
        fig.savefig(fname, dpi=150, bbox_inches='tight')
        # plt.show()


def add_event(model_string, event_str):
    model_string = model_string[:-5]
    event_str = "    " + event_str
    model_string += event_str
    model_string += "\nend;"
    return model_string


def add_serum_starve_event_remove_growth_factors(model_string):
    """
    Set growth factors to 0 at t=70.25 minutes
    :param model_string:
    :return:
    """
    event_str = "SerumStarveRemoveGrowthFactors: at (time>70.25): GrowthFactors=0;\n"
    model_string = add_event(model_string, event_str)
    return model_string

def add_serum_starve_event_remove_basal_TGFb(model_string):
    """
    Set growth factors to 0 at t=70.25 minutes
    :param model_string:
    :return:
    """
    event_str = "SerumStarveRemoveTGFb: at (time>70.25): TGFb=0;\n"
    model_string = add_event(model_string, event_str)
    return model_string


def add_TGFb_event(model_string):
    """
    Put 1 TGFb into simulation at t=4275
    :param model_string:
    :return:
    """
    event_str = "AddTGFb: at (time>71.250): TGFb=1;\n"
    model_string = add_event(model_string, event_str)
    return model_string


def add_MK_event(model_string, time):
    """
    Put 1 TGFb into simulation at t=4275
    :param model_string:
    :return:
    """
    event_str = "AddMK2206: at (time>{}): MK2206=1;\n".format(time)
    model_string = add_event(model_string, event_str)
    return model_string


def add_AZD_event(model_string, time):
    """
    Put 1 TGFb into simulation at t=4275
    :param model_string:
    :return:
    """
    event_str = "AddAZD: at (time>{}): AZD=1;\n".format(time)
    model_string = add_event(model_string, event_str)
    return model_string


def make_condition(model_string, condition):
    """
    Take the model string and create condition

    :param model_string: antimony model string
    :param GF: starting amount of GrowthFactors
    :param TGF: starting amount of tgf
    :param pretreatment: either 'AZD', or 'MK2206'. This is added as an event at the time specified by pretreatment time
    :param pretreatment_time: The time at which the variable specified by 'pretreatment' is added
    :param EV: Starting amount of Everolimus
    :param serum_starve_event: Boolean, whether to remove serum, aka GrowthFactors by event
    :param TGFb_event: Boolean. Whether to add 1 unit TGF at t=-45min (71.75h)
    :param open_with_copasi: open the model with copasi
    :return:
    """
    if condition not in AZD_CONDITIONS.keys() + MK_CONDITIONS.keys():
        raise ValueError('No key. These keys "{}"'.format(AZD_CONDITIONS.keys()+MK_CONDITIONS.keys()))
    try:
        GF, pretreatment, pretreatment_time, EV, serum_starve_event, TGFb_event = AZD_CONDITIONS[condition]
    except KeyError:
        GF, pretreatment, pretreatment_time, EV, serum_starve_event, TGFb_event = MK_CONDITIONS[condition]

    if serum_starve_event:
        model_string = add_serum_starve_event_remove_basal_TGFb(model_string)
        model_string = add_serum_starve_event_remove_growth_factors(model_string)

    if TGFb_event:
        model_string = add_TGFb_event(model_string)
    if pretreatment is not None:
        if pretreatment_time != 0:
            if pretreatment == 'AZD':
                model_string = add_AZD_event(model_string, pretreatment_time)
            elif pretreatment == 'MK2206':
                model_string = add_MK_event(model_string, pretreatment_time)
            else:
                raise ValueError

    mod = te.loada(model_string)
    # mod.model['init([TGFb])'] = TGF
    mod.Everolimus = EV
    mod.GrowthFactors = GF

    if pretreatment_time == 0:
        if pretreatment == 'AZD':
            mod.AZD = 1
        elif pretreatment == 'MK2206':
            mod.MK2206 = 1

    return mod.getCurrentAntimony()


def simulate_condition(model_string, condition):
    """
    :param model_string: antimony model string
    :param GF: starting amount of GrowthFactors
    :param pretreatment: either 'AZD', or 'MK2206'. This is added as an event at the time specified by pretreatment time
    :param pretreatment_time: The time at which the variable specified by 'pretreatment' is added
    :param EV: Starting amount of Everolimus
    :param serum_starve_event: Boolean, whether to remove serum, aka GrowthFactors by event
    :param TGFb_event: Boolean. Whether to add 1 unit TGF at t=-45min (71.75h)
    :return:
    """
    if condition not in AZD_CONDITIONS.keys() + MK_CONDITIONS.keys():
        raise ValueError

    model_string = make_condition(model_string, condition)
    mod = te.loada(model_string)
    ## Add 1 to intervals for 0 indexed python
    start = 0
    end = 75
    num = 751
    time = numpy.linspace(start, end, num)

    if 72.0 not in time:
        raise ValueError("Simulating conditions requires an output at exactly "
                         "72h after the first inhibitor. The current time steps "
                         "do not 'land' on 72.0 exactly - so pick different "
                         "simulation parameters. Also, remember that Python is "
                         "0 indexed, so try adding 1 to the number of "
                         " data points collected. To help, here is the current"
                         " time vector \n\n{}".format(time))

    res = mod.simulate(start, end, num, ['time'] + MODEL_SPECIES)

    df = pandas.DataFrame(res, columns=['time'] + MODEL_SPECIES)
    return df


def simulate_conditions(model_str=None, type='azd'):
    """
    Takes output from simulation and plot
    :param df:
    :return:
    """
    if model_str is None:
        model_str = cross_talk_model_antstr()

    seaborn.set_context('talk', font_scale=2)
    if type not in ['AZD', 'MK2206']:
        raise ValueError

    if type == 'AZD':
        conditions = AZD_CONDITIONS
    elif type == 'MK2206':
        conditions = MK_CONDITIONS
    else:
        raise ValueError

    # from multiprocessing import Process

    dct = OrderedDict()
    for k, v in conditions.items():
        # p = Process(target=simulate_condition, args=(tuple([cross_talk_model_antstr()] + conditions[k])))
        # p.start()
        # p.join()
        df = simulate_condition(model_str, k)
        df = df[df['time'] == 72.0]
        if df.empty:
            raise ValueError("Condition '{}' with values '{}' produces an "
                             "empty data frame".format(k, v))

        dct[k] = df

    df = pandas.concat(dct)
    df.index = df.index.droplevel(1)
    df = df.drop('time', axis=1)


    return df


def simulate_conditions_and_plot_as_bargraph(y, type='AZD'):
    """
    Takes output from simulation and plot
    :param df:
    :return:
    """

    df = simulate_conditions(type=type)
    if y not in df.columns:
        raise ValueError('The variable "{}" does not appear in your model'.format(y))
    df = df.reset_index()
    df = df.rename(columns={'index': 'Condition'})

    fig = plt.figure()
    seaborn.barplot(x='Condition', y=y, data=df,
                    palette=['yellow'] * 2 + ['white'] * 1 + ['red'] * 4 + ['green'] * 4,
                    edgecolor='black', linewidth=2,
                    order=MK_CONDITIONS.keys() if type == 'MK2206' else AZD_CONDITIONS.keys()
                    )
    plt.title(y)
    seaborn.despine(fig, top=True, right=True)
    plt.xticks(rotation=90, fontsize=18)
    plt.yticks(fontsize=18)

    if type == 'AZD':
        dire = os.path.join(graphs_directory, 'AZDSimulations')
    else:
        dire = os.path.join(graphs_directory, 'MKSimulations')
    os.makedirs(dire) if not os.path.isdir(dire) else None
    fname = os.path.join(dire, "{}.png".format(y))
    fig.savefig(fname, dpi=300, bbox_inches='tight')
    print('line 485: Figure saved to "{}"'.format(fname))



def simulate_all_conditions_and_plot_as_bargraphs():
    """
    Plot all bar graphs from simulations
    :return:
    """
    model_species = MODEL_SPECIES
    for i in model_species:
        simulate_conditions_and_plot_as_bargraph(i, 'AZD')
        simulate_conditions_and_plot_as_bargraph(i, 'MK2206')


def simulate_timecourse(type='MK2206'):
    """
    Events are used throughout the simulations.
    This function is for visualising the model inputs
    to visualise whether the events are correctly configured.
    :return:
    """
    if type not in ['AZD', 'MK2206']:
        raise ValueError

    if type == 'AZD':
        conditions = AZD_CONDITIONS
    elif type == 'MK2206':
        conditions = MK_CONDITIONS
    else:
        raise ValueError

    dct = OrderedDict()
    for k, v in conditions.items():
        df = simulate_condition(cross_talk_model_antstr(), k)
        dct[k] = df
    # print(dct)

    df = pandas.concat(dct)
    df.index = df.index.droplevel(1)
    df = df.set_index('time', append=True)
    df.index = df.index.rename(['condition', 'time'])

    return df


def plot_timecourse_single(vars, condition=None):
    if condition not in AZD_CONDITIONS.keys() + MK_CONDITIONS.keys():
        raise ValueError

    df = simulate_timecourse()
    for i in vars:
        if i not in df.columns:
            raise ValueError("Variable '{}' is not in '{}'".format(i, df.columns))
    df = df[vars]

    df = df.reset_index()

    if condition is None:
        condition = df['condition'].unique()

    df = df[df['condition'] == condition]

    fig_dct = {}

    for v in range(len(vars)):
        fig = plt.figure()
        seaborn.despine(fig=fig, top=True, right=True)
        plt.title("{} in {} condition".format(vars[v], condition))
        plt.xlabel('Time(h)')
        plt.ylabel('nmol L$^{-1}$')
        plt.plot(df['time'], df[vars[v]], label=vars[v])

        fig_dct["{}_{}".format(vars[v], condition)] = fig

    return fig


def plot_timecourse_multiplot(vars, condition=None, multifig_shape=None, AZD_or_MK='AZD'):


    # if condition not in AZD_CONDITIONS.keys() + MK_CONDITIONS.keys():
    #     raise ValueError

    if multifig_shape is not None:
        assert len(multifig_shape) == 2
        assert isinstance(multifig_shape, (tuple, list))

    if multifig_shape is None:
        multifig_shape = [vars, 1]

    if AZD_or_MK not in ['AZD', 'MK2206']:
        raise ValueError

    if AZD_or_MK == 'MK2206':
        if condition not in MK_CONDITIONS.keys():
            raise ValueError

    elif AZD_or_MK == 'AZD':
        if condition not in AZD_CONDITIONS.keys():
            raise ValueError('"{}" condition is not in "{}"'.format(condition, AZD_CONDITIONS.keys()))

    df = simulate_timecourse(AZD_or_MK)


    for i in vars:
        if i not in df.columns:
            raise ValueError("Variable '{}' is not in '{}'".format(i, df.columns))

    df = df[vars]

    df = df.reset_index()

    if condition is None:
        condition = df['condition'].unique()

    df = df[df['condition'] == condition]

    fig = plt.figure()
    ax = plt.subplot2grid(multifig_shape, [0, 0])


    cols = seaborn.color_palette("hls", len(vars))
    cols = iter(cols)

    for v in range(len(vars)):
        seaborn.set_context('talk')
        sub_ax = plt.subplot(*multifig_shape + [v+1])
        plt.setp(sub_ax.get_xticklabels(), visible=False)
        # box = sub_ax.get_position()
        # sub_ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

        seaborn.despine(ax=sub_ax, top=True, right=True)
        plt.yticks(fontsize=12)
        plt.suptitle(condition)
        plt.plot(df['time'], df[vars[v]], label=vars[v], color=next(cols))
        plt.legend(loc=(0.75, 0.5), bbox_to_anchor=(1.0, 0.5), fontsize=10)
        plt.xticks(numpy.arange(0, 72, 5.0), fontsize=12, rotation=90)

        plt.subplots_adjust(bottom=0.2, left=0.2, right=0.7)
        if v == len(vars)-1:
            plt.xlabel('Time(h)')

    plt.setp(sub_ax.get_xticklabels(), visible=True)

    fig.text(0.05, 0.5, 'nmol L$^{-1}$', ha='center', va='center', rotation='vertical')

    return fig


def simulate_inputs_only(AZD_or_MK='AZD'):
    if AZD_or_MK not in ['AZD', 'MK2206']:
        raise ValueError

    dire = os.path.join(graphs_directory, 'InputSimulations')
    os.makedirs(dire) if not os.path.isdir(dire) else None

    if AZD_or_MK == 'AZD':
        conditions = AZD_CONDITIONS.keys()
    elif AZD_or_MK == 'MK2206':
        conditions = MK_CONDITIONS.keys()
    else:
        raise ValueError

    conditions = list(set(conditions))

    for cond in conditions:
        fig = plot_timecourse_multiplot(['TGFb', 'MK2206', 'Everolimus', 'AZD', 'GrowthFactors'],
                                        condition=cond, multifig_shape=[5, 1], AZD_or_MK=AZD_or_MK)

        fname = os.path.join(dire, "{}.png".format(cond))
        fig.savefig(fname, dpi=300, bobx_inches='tight')
        print("Saved to '{}'".format(fname))




def open_condition_with_copasi(model_string, condition):
    """
    Create condition with model string and open with copasi,
    saving to folder while your at it
    :param condition:
    :return:
    """
    mod_string = make_condition(model_string, condition)
    copasi_dir = os.path.join(working_directory, 'copasi_models')
    if not os.path.isdir(copasi_dir):
        os.makedirs(copasi_dir)
    fname = os.path.join(copasi_dir, "{}.cps".format(condition))
    copasi_mod = load_model_with_pyco(mod_string, fname)
    tasks.TimeCourse(copasi_mod, start=0, end=72, step_size=0.1, intervals=720, run=False)
    copasi_mod.open()
    return copasi_mod


def configure_parameter_estimation(model_string, condition):
    mod_string = make_condition(model_string, condition)
    copasi_dir = os.path.join(working_directory, 'copasi_models')
    if not os.path.isdir(copasi_dir):
        os.makedirs(copasi_dir)
    fname = os.path.join(copasi_dir, "{}.cps".format(condition))
    copasi_mod = load_model_with_pyco(mod_string, fname)
    tasks.TimeCourse(copasi_mod, start=0, end=72, step_size=0.1, intervals=720, run=False)
    f1 = r'D:\MesiSTRAT\CrossTalkModel\copasi_models\fakeSteadyStateData.txt'
    f2 = r'D:\MesiSTRAT\CrossTalkModel\copasi_models\Smad2.txt'
    assert os.path.isfile(f1)
    assert os.path.isfile(f2)
    # copasi_mod.open()
    PE = tasks.ParameterEstimation(copasi_mod, [f1, f2], metabolites=[], run_mode=False,
                                   experiment_type=['steadystate', 'steadystate'],
                                   overwrite_config_file=False)
    PE.write_config_file()
    PE.setup()
    return PE.model




def simulate_model_component_timecourse(vars, cond, filename=None, **kwargs):
    """
    Simulate variables in vars in conditions in cond and plot on one plot
    :param var:
    :param cond:
    :return:
    """
    seaborn.set_context('talk', font_scale=2)
    if isinstance(vars, str):
        vars = [vars]

    if isinstance(cond, str):
        cond = [cond]

    for i in vars:
        if i not in MODEL_SPECIES:
            raise ValueError("model species '{}' is not in '{}'".format(i, MODEL_SPECIES))

    if kwargs.get('figsize') is not None:
        figsize = kwargs['figsize']
        del kwargs['figsize']
    else:
        figsize = (2*len(cond), 3*len(cond))

    cols = seaborn.color_palette("hls", len(vars))*len(cond)
    cols = iter(cols)

    fig = plt.figure(figsize=figsize)
    for i, c in enumerate(cond):
        if c not in list(MK_CONDITIONS.keys()) + list(AZD_CONDITIONS.keys()):
            raise ValueError('condition "{}" not in "{}"'.format(c, MK_CONDITIONS.keys() + AZD_CONDITIONS.keys()))

        model_str = make_condition(cross_talk_model_antstr(), c)
        mod = te.loada(model_str)
        res = mod.simulate(0, 72, 73, ['time'] + vars)
        gs = GridSpec(len(cond), 1, wspace=0.3)
        ax = []
        for j, v in enumerate(vars):
            ax.append(fig.add_subplot(gs[i, 0]))
            ax[-1].plot(res['time'], res[v], label=v, color=next(cols), linewidth=6, **kwargs)

        plt.axvline(0.0  ,   linestyle='--', linewidth=2, color='black',  alpha=0.4)
        plt.axvline(24.0 ,   linestyle='--', linewidth=2, color='black',  alpha=0.4)
        plt.axvline(48.0 ,   linestyle='--', linewidth=2, color='black',  alpha=0.4)
        plt.axvline(70.25,   linestyle='--', linewidth=2, color='green',  alpha=0.4)
        plt.axvline(71.25,   linestyle='--', linewidth=2, color='purple', alpha=0.4)
        plt.axvline(72.0 ,   linestyle='--', linewidth=2, color='black',  alpha=0.4)

        plt.setp(ax[-1].get_xticklabels(), visible=False)
        seaborn.despine(ax=ax[-1], top=True, right=True)
        plt.title(c)
    if i == len(cond):
        plt.xlabel('Time(h)')

    fig.text(0.025, 0.5, 'nmol L$^{-1}$', ha='center', va='center', rotation='vertical')
    plt.legend(loc=[1, 1])
    plt.setp(ax[-1].get_xticklabels(), visible=True)

    dire = os.path.join(graphs_directory, 'TimeSeriesSimulations')
    os.makedirs(dire) if not os.path.isdir(dire) else None
    if filename is None:
        fname = os.path.join(dire, reduce(lambda x, y: "{}__{}".format(x, y), cond) + '.png')
    else:
        fname = os.path.join(dire, filename)

    plt.savefig(fname, dpi=250, bbox_inches='tight')


def get_parameters_from_copasi_in_antimony_format(condition):
    copasi_dir = os.path.join(working_directory, 'copasi_models')
    fname = os.path.join(copasi_dir, "{}.cps".format(condition))
    sbml = os.path.join(copasi_dir, "{}.sbml".format(condition))
    pm = model.Model(fname)
    pm.to_sbml(sbml)
    mod = te.loadSBMLModel(sbml)
    print(te.getCurrentAntimony(mod))

def get_model_parameters(mod):
    return dict(zip(mod.getGlobalParameterIds(), mod.getGlobalParameterValues()))

class Inequality(object):
    def __init__(self, left, operator, right, name):
        """

        :param left: [cond, species]
        :param operator: '>' or '<'
        :param right: [cond, species]
        """
        self.left = left
        self.operator = operator
        self.right = right
        self.name = name

    def evaluate(self, df):
        left = df.loc[self.left[0], self.left[1]]
        right = df.loc[self.right[0], self.right[1]]
        if self.operator == '>':
            return left > right
        elif self.operator == '<':
            return left < right

    def __str__(self):
        return 'Inequality(left="{}", operator="{}", right="{}", name="{}")'.format(
            self.left, self.operator, self.right, self.name
        )

    def __repr__(self):
        return self.__str__()


class InequalityGroup(object):
    def __init__(self, inequalities):
        self.inequalities = inequalities

        for ineq in self.inequalities:
            self.__dict__[ineq.name] = ineq

    def __iter__(self):
        return self

    def __getitem__(self, item):
        return self.__dict__[item]

    def __setitem__(self, key, value):
        self.__dict__[key] = value

    def __delitem__(self, key):
        del self.__dict__[key]

    def next(self):
        count = 0
        if count > self.__len__():
            raise StopIteration
        else:
            count += 1
            return self.inequalities.keys()[count]

    def inequality_keys(self):
        return cycle(range(len(self)))

    @property
    def names(self):
        return [i.name for i in self.inequalities]

    def evaluate(self, df):
        eval_dct = OrderedDict()
        for ineq in self.inequalities:
            if not isinstance(ineq, Inequality):
                raise ValueError

            ev = ineq.evaluate(df)
            eval_dct[ineq.name] = ev

        return eval_dct

    def __len__(self):
        return len(self.inequalities)


class OptimizeQualitative(object):
    pandas.set_option('precision', 24)
    getcontext().prec = 24

    @property
    def condition_list(self):
        return cycle(self.inequality_group.names)

    def __init__(self, model_string, free_parameters,
                 inequality_group, exp_data, delta, iterations=10
                 ):
        self.model_string = model_string
        self.free_parameters = free_parameters
        self.exp_data = exp_data
        self.inequality_group = inequality_group
        self.delta = delta
        self.fname_param = os.path.join(working_directory, 'parameters.csv')
        self.fname_prob = os.path.join(working_directory, 'prob.csv')
        self.iterations = iterations

        self.mod = self.load_model()

    def load_model(self):
        return te.loada(self.model_string)

    @staticmethod
    def simulate_condition(model_string, condition):
        if condition not in AZD_CONDITIONS.keys() + MK_CONDITIONS.keys():
            raise ValueError

        model_string = make_condition(model_string, condition)
        mod = te.loada(model_string)
        ## Add 1 to intervals for 0 indexed python
        start = 0
        end = 75
        num = 751
        time = numpy.linspace(start, end, num)

        if 72.0 not in time:
            raise ValueError("Simulating conditions requires an output at exactly "
                             "72h after the first inhibitor. The current time steps "
                             "do not 'land' on 72.0 exactly - so pick different "
                             "simulation parameters. Also, remember that Python is "
                             "0 indexed, so try adding 1 to the number of "
                             " data points collected. To help, here is the current"
                             " time vector \n\n{}".format(time))

        res = mod.simulate(start, end, num, ['time'] + MODEL_SPECIES)

        df = pandas.DataFrame(res, columns=['time'] + MODEL_SPECIES)
        df = df[df['time'] == 72.0]
        if df.empty:
            raise ValueError("Condition '{}' produces an "
                             "empty data frame".format(condition))
        return df

    def _simulate_conditions(self, model_str, conditions=None):
        """
        Takes output from simulation and plot
        :param df:
        :return:
        """
        if conditions is None:
            conditions = list(set(AZD_CONDITIONS.keys() + MK_CONDITIONS.keys()))

        dct = OrderedDict()
        for cond in conditions:
            df = self.simulate_condition(model_str, cond)
            df = df[df['time'] == 72.0]
            dct[cond] = df

        df = pandas.concat(dct)
        df.index = df.index.droplevel(1)
        df = df.drop('time', axis=1)
        return df

    def _evaluate_inequalities(self, df):
        self.inequality_group.evaluate(df)
        return self.inequality_group.memory

    def _check_if_all_conditions_eval_to_true(self, inequality_dct):
        if inequality_dct is None:
            return None

        if not isinstance(inequality_dct, (dict, None)):
            raise TypeError("Argument should be a dictionary. Got a '{}' instead.".format(type(inequality_dct)))

        if all(inequality_dct.values()):
            ## return if all values are True
            print('all conditions are met')
            print(self.free_parameters)
            return True
        else:
            return False

    def _peterb_parameter(self, mod, param):
        old = self.free_parameters[param]
        peterb_amount = numpy.random.normal(0, 0.1 * old, 1)[0]
        assert peterb_amount < old, "peterb amount must be smaller than old value otherwise risk getting negative values"
        new_parameter_value = old + peterb_amount
        assert new_parameter_value > 0, "Be positive"
        setattr(mod, param, new_parameter_value)
        print('old parameter value is "{}". New is "{}"'.format(old, new_parameter_value))
        return mod, new_parameter_value

    def _load_parameters_from_dict(self, free_params):
        for i in free_params:
            setattr(self.mod, i, free_params[i])

    def fit(self):
        """

        :return:
        """
        ## if all conditions are satisified, return
        df = self._simulate_conditions(self.model_string)
        ineq_memory = self._evaluate_inequalities(df)
        if self._check_if_all_conditions_eval_to_true(ineq_memory['new']):
            self.save_parameters_and_prob_to_file()
            return self.free_parameters, self.probability_matrix

        # mod = te.loada(self.model_string)
        for i in range(self.iterations):
            print('\nIteration {}'.format(i))
            ## load self.free_parameters back into the model
            self._load_parameters_from_dict(self.free_parameters)
            # print(zip(self.mod.getGlobalParameterIds(), self.mod.getGlobalParameterValues()))
            current_condition = next(self.condition_list)
            ## choose one parameter to peterb, based on probabilities
            choice = numpy.random.choice(self.free_parameters.keys(),
                                         p=self.probability_matrix[current_condition].values)

            # print('Current choice is "{}", value "{}"'.format(choice, getattr(self.mod, choice)))
            self.mod, new_parameter_value = self._peterb_parameter(self.mod, choice)
            df = self._simulate_conditions(self.mod.getCurrentAntimony())
            ineq_memory = self._evaluate_inequalities(df)

            if self._check_if_all_conditions_eval_to_true(ineq_memory['new']):
                self.save_parameters_and_prob_to_file()
                return self.free_parameters, self.probability_matrix

            for k in ineq_memory['old'].keys():
                if ineq_memory['old'][k] != ineq_memory['new'][k]:
                    print ('ineq "{}" changed from "{}" to "{}"'.format(k, ineq_memory['old'][k], ineq_memory['new'][k]))

            print('old count: {}'.format(Counter(ineq_memory['old'].values())))
            print('new count: {}'.format(Counter(ineq_memory['new'].values())))

            old_parameter_value = self.free_parameters[choice]

            for cond, boolean in ineq_memory['new'].items():
                # print ('cond={}, cond {}, ineq_memory old "{}", ineq_memory new "{}" eq {}'.format(
                #     cond,
                #     boolean,
                #     ineq_memory['old'][cond],
                #     ineq_memory['new'][cond],
                #     (ineq_memory['old'][cond] == False) and (ineq_memory['old'][cond] == True)
                # ))
                if (ineq_memory['old'][cond] == False) and (ineq_memory['new'][cond] == True):
                    print("Positively reinforcing '{}' parameter as changing it from '{}' to '{}' "
                          "changed the '{}' condition from False to True".format(
                        choice, old_parameter_value, new_parameter_value, cond
                    ))
                    self._update_probabilities(cond, choice, 'positive_reinforcement')
                    self.free_parameters[choice] = new_parameter_value

                elif (ineq_memory['old'][cond] == True) and (ineq_memory['new'][cond] == False):
                    print("Negatively reinforcing '{}' parameter as changing it from '{}' to '{}' "
                          "changed the '{}' condition from True to False".format(
                        choice, old_parameter_value, new_parameter_value, cond
                    ))
                    self._update_probabilities(cond, choice, 'negative_reinforcement')

        self.save_parameters_and_prob_to_file()
        print(self.free_parameters)
        return self.free_parameters, self.probability_matrix


class Mapping(object):
    def __init__(self, experimental, simulation):
        self.experimental = experimental
        self.simulation = simulation

    def __str__(self):
        return ("({}, {})".format(self.experimental, self.simulation))


# class MappingGroup(object):
#     def __init__(self, mappings):
#         self.mappings = mappings
#
#     def __iter__(self):
#         return self
#
#     def __nex


class RSS(object):
    def __init__(self, simulated_data, experimental_data, mappings):
        self.simulated_data = simulated_data
        self.experimental_data = experimental_data
        self.mappings = mappings

        # self.residuals = self._compute()
        # self.rss = self.residuals.sum().sum()

    def _compute(self):
        dct = OrderedDict()
        for i in self.mappings:
            exp = i.experimental
            sim = i.simulation
            # print(sim, exp)
            experimental = self.experimental_data[exp]
            # print(experimental)
            simulated = self.simulated_data[sim]
            # print(simulated)
            dct[sim] = (experimental - simulated)**2
        return pandas.concat(dct, axis=1)


class QualitativeObj(object):
    """

    """
    def __init__(self, simulated_data, inequalities):
        self.simulated_data = simulated_data
        self.inequalities = inequalities
        # print(self.simulated_data)
        # print(self.inequalities)

        self.residuals = self._compute()
        self.residual_values = numpy.array(self.residuals.values())
        self.negative_loss = sum([i for i in self.residual_values if i < 0])
        self.positive_loss = sum([i for i in self.residual_values if i > 0])

    def __str__(self):
        return 'QualitativeObj(positive_loss={}, negative_loss={})'.format(self.positive_loss, self.negative_loss)

    def __repr__(self):
        return self.__str__()

    def _compute(self):
        """
        Value of 0 means that inequality is satisified
        +ve value means that left is bigger than right
        -ve value means that left is smaller than right
        :return:
        """
        obj = OrderedDict()
        dct = self.inequalities.evaluate(self.simulated_data)
        for k, v in dct.items():

            if v == True:
                obj[k] = 0
            else:
                ineq = self.inequalities[k]
                ## left is always small
                if ineq.operator == '<':
                    left = self.simulated_data.loc[ineq.right[0], ineq.right[1]]
                    right = self.simulated_data.loc[ineq.left[0], ineq.left[1]]
                else:
                    left = self.simulated_data.loc[ineq.left[0], ineq.left[1]]
                    right = self.simulated_data.loc[ineq.right[0], ineq.right[1]]

                obj[k] = right - left
        return obj



class RandomSimulation(object):
    """

    I am trying to produce some simulation data for anaysis with a neural net.

    To do this I need to:
        1) simulate data with parameter set
        2) evaluate truth of inequalities
            - Also compute the qualitative objecctive function
            - 0 if True and +/- if not
        3) return dataset including the parameters and the inequality results.


    """
    def __init__(self, model_string, exp_data, inequalities, mappings, free_parameters, iterations):
        self.model_string = model_string
        self.exp_data = exp_data
        self.inequalities = inequalities
        self.mappings = mappings
        self.free_parameters = free_parameters
        self.iterations = iterations

        self.model = self._load_model()

        self.params = {}
        self.resid = {}

        self.fname_params = os.path.join(working_directory, 'parameters.csv')
        self.fname_resid = os.path.join(working_directory, 'residuals.csv')

    def _load_model(self):
        return te.loada(self.model_string)

    def _simulate_condition(self, condition):
        if condition not in AZD_CONDITIONS.keys() + MK_CONDITIONS.keys():
            raise ValueError

        model_string = make_condition(self.model.getCurrentAntimony(), condition)
        mod = te.loada(model_string)
        ## Add 1 to intervals for 0 indexed python
        start = 0
        end = 72
        num = 73
        time = numpy.linspace(start, end, num)

        if 72.0 not in time:
            raise ValueError("Simulating conditions requires an output at exactly "
                             "72h after the first inhibitor. The current time steps "
                             "do not 'land' on 72.0 exactly - so pick different "
                             "simulation parameters. Also, remember that Python is "
                             "0 indexed, so try adding 1 to the number of "
                             " data points collected. To help, here is the current"
                             " time vector \n\n{}".format(time))

        res = mod.simulate(start, end, num, ['time'] + MODEL_SPECIES)

        df = pandas.DataFrame(res, columns=['time'] + MODEL_SPECIES)
        df = df[df['time'] == 72.0]
        if df.empty:
            raise ValueError("Condition '{}' produces an "
                             "empty data frame".format(condition))
        return df

    def _simulate_conditions(self, conditions=None):
        """
        Takes output from simulation and plot
        :param df:
        :return:
        """
        if conditions is None:
            conditions = list(set(AZD_CONDITIONS.keys() + MK_CONDITIONS.keys()))

        dct = OrderedDict()
        for cond in conditions:
            df = self._simulate_condition(cond)
            df = df[df['time'] == 72.0]
            dct[cond] = df

        df = pandas.concat(dct)
        df.index = df.index.droplevel(1)
        df = df.drop('time', axis=1)
        return df

    def _compute_obj_fun(self, sim_data):
        sim_names = [i.simulation for i in self.mappings]
        exp_names = [i.experimental for i in self.mappings]
        sim_data = sim_data[sim_names]
        qual = QualitativeObj(sim_data, self.inequalities)
        return qual

    @staticmethod
    def _load_parameters_from_dict(mod, free_params):
        if not isinstance(free_params, dict):
            raise TypeError
        for i in free_params:
            setattr(mod, i, free_params[i])
        return mod

    def fit1(self, i):
        # params = {}
        # resid = {}
        new_params = numpy.random.uniform(0.01, 1000, len(self.free_parameters))
        new_params = OrderedDict(zip(self.free_parameters.keys(), new_params))
        self.model = self._load_parameters_from_dict(self.model, new_params)
        sim_data = self._simulate_conditions()
        obj = self._compute_obj_fun(sim_data)
        # params[i] = new_params
        # resid[i] = obj.residuals
        return new_params, obj.residuals#params, resid

    def to_file(self, params, resid):
        if os.path.isfile(self.fname_params):
            with open(self.fname_params, 'a') as f:
                params.to_csv(f, header=False)
        else:
            params.to_csv(self.fname_params)

        if os.path.isfile(self.fname_resid):
            with open(self.fname_resid, 'a') as f:
                resid.to_csv(f, header=False)
        else:
            resid.to_csv(self.fname_resid)

        print('parameters saved to "{}"'.format(self.fname_params))
        print('residuals saved to "{}"'.format(self.fname_resid))


    def fit(self):
        import time
        resid = {}
        params = {}

        try:
            for i in (j for j in range(self.iterations)):
                print('iteration {}'.format(i))
                now = time.clock()
                params[i], resid[i] = self.fit1(i)
                print("loop took '{}' seconds".format(time.clock() - now))

        except KeyboardInterrupt:
            print('Warning: Keyboard Interupt. Saving progress.')

        params = pandas.DataFrame(params).transpose()
        resid = pandas.DataFrame(resid).transpose()

        self.to_file(params, resid)

        return params, resid





if __name__ == '__main__':

    working_directory = r'D:\MesiSTRAT\CrossTalkModel'
    copasi_filename = os.path.join(working_directory, 'KatrinesTopology.cps')
    graphs_directory = os.path.join(working_directory, 'SimulationGraphs')

    """
    Set flags to determine which part of the script will run
    """

    def set_flags():
        """
        Dummy function to get clickable button on
        PyCharm structure tab to take me directly to
        the flags
        :return:
        """

    if GET_PARAMETERS_FROM_COPASI:
        get_parameters_from_copasi_in_antimony_format(MODEL_CODE)

    if OPEN_CONDITION_WITH_COPASI:
        open_condition_with_copasi(cross_talk_model_antstr(), MODEL_CODE)

    if CONFIGURE_PARAMETER_ESTIMATION:
        m = configure_parameter_estimation(cross_talk_model_antstr(), MODEL_CODE)
        m.open()

    phos = ['pErk', 'pAkt', 'pSmad2', 'pRaf', 'ppMek', 'ppErk',
            'pPI3K', 'pPI3K', 'pmTORC1', 'pS6K']
    erk = ['Erk', 'pErk', 'ppErk']

    # pSmad2  =   ['pSmad2Tot']

    if SIMULATE_TIME_SERIES:
        for i in CURRENT_SPECIES:
            simulate_model_component_timecourse([i], AZD_CONDITIONS.keys(), filename='AZD_'+i)
            simulate_model_component_timecourse([i], MK_CONDITIONS.keys(), filename='MK_'+i)

    if SIMULATE_BAR_GRAPHS:
        # for i in ['AZD', 'MK2206']:
        for i in ['MK2206', 'AZD']:
            for j in CURRENT_SPECIES:
                simulate_conditions_and_plot_as_bargraph(j, i)



    if SIMULATE_INPUTS:
        simulate_inputs_only('AZD')
        simulate_inputs_only('MK2206')

    if not os.path.isdir(graphs_directory):
        os.makedirs(graphs_directory)

    if GET_ODES_WITH_ANTIMONY:
        mod = te.loada(cross_talk_model_antstr())
        odes = te.utils.misc.getODEsFromModel(mod)
        print(odes)

    if GET_MODEL_AS_SBML:
        sbml = te.antimonyToSBML(cross_talk_model_antstr())
        print(sbml)

    if DOSE_RESPONSE_TGFB:
        fig = dose_response(cross_talk_model_antstr(), 'TGFb', 0.01, 1000, 100, ['TGFbR'])
        plt.show()


    if QUALITATIVE_FITTING:
        model_string = cross_talk_model_antstr()

        azd_data = os.path.join(working_directory, r'data\HardCopy\AZD_calculations - v3.xlsx')
        mk_data = os.path.join(working_directory, r'data\HardCopy\MK2206_calculations - v3.xlsx')

        assert os.path.isfile(azd_data)
        assert os.path.isfile(mk_data)

        azd_data = pandas.read_excel(azd_data, sheet_name='AZDAverage', index_col='Condition')
        mk_data = pandas.read_excel(mk_data, sheet_name='MKAverage', index_col='Condition')

        mk_data = mk_data.iloc[3:]

        exp_data = pandas.concat([azd_data, mk_data])

        free_parameters = OrderedDict({
            'kmTORC1Phos_ki': 0.001,
            'kPI3KPhosByTGFbR_kcat': 50.0,
            'kAktDephos_Vmax': 31.1252344504785,
            'kPI3KDephosByErk': 5.014,
            'kS6KPhosBymTORC1_kcat': 2.77975221288272,
            'kPI3KPhosByGF': 0.239474698704283,
            'kPI3KDephosByS6K': 25.0,
            'kErkPhos_kcat1': 85.0103161451182,
            'kmTORC1Dephos_Vmax': 1.0,
            'kS6KDephos_Vmax': 50.0,
            'kAktPhos_kcat': 2.9215,
            'kRafPhos_ki': 3.5,
            'kRafPhosByTGFbR_kcat': 265.0,
            'kRafPhosByPI3K_kcat': 50.0,
            'kMekPhos_kcat1': 149.5209856,
            'kMekPhos_ki1': 0.25,
            'kTGFbOn': 0.100647860357268,
            'kSmad2PhosByAkt_kcat': 1.0,
            'kpSmad2Dephos_Vmax': 58.8712661228653,
            'kAktPhos_ki': 0.01,
            'kmTORC1Phos_kcat': 0.1,
        })
        # for i in free_parameters:
        #     free_parameters[i] = numpy.random.uniform(0, 100, 1)[0]

        delta = 0.01
        akt_ineq1  = Inequality(['E',        'pAkt'],  '>', ['D',        'pAkt'], 'pAkt_1')
        akt_ineq2  = Inequality(['E_A_72',   'pAkt'],  '>', ['E',        'pAkt'], 'pAkt_2')
        akt_ineq3  = Inequality(['A_72',     'pAkt'],  '<', ['E',        'pAkt'], 'pAkt_3')
        akt_ineq4  = Inequality(['A_72',     'pAkt'],  '<', ['E_A_72',   'pAkt'], 'pAkt_4')
        akt_ineq5  = Inequality(['E_M_72',   'pAkt'],  '<', ['D',        'pAkt'], 'pAkt_5')
        akt_ineq6  = Inequality(['E_M_72',   'pAkt'],  '<', ['T',        'pAkt'], 'pAkt_6')
        akt_ineq7  = Inequality(['E_M_72',   'pAkt'],  '<', ['E',        'pAkt'], 'pAkt_7')
        akt_ineq8  = Inequality(['M_72',     'pAkt'],  '<', ['D',        'pAkt'], 'pAkt_8')
        akt_ineq9  = Inequality(['M_72',     'pAkt'],  '<', ['T',        'pAkt'], 'pAkt_9')
        akt_ineq10 = Inequality(['M_72',     'pAkt'],  '<', ['E',        'pAkt'], 'pAkt_10')

        akt = [
            akt_ineq1,
            akt_ineq2,
            akt_ineq3,
            akt_ineq4,
            akt_ineq5,
            akt_ineq6,
            akt_ineq7,
            akt_ineq8,
            akt_ineq9,
            akt_ineq10,
        ]

        erk_ineq1 = Inequality(['E',        'ppErk'],  '>', ['D',        'ppErk'], 'ppErk_1')
        erk_ineq2 = Inequality(['E',        'ppErk'],  '>', ['T',        'ppErk'], 'ppErk_2')
        erk_ineq3 = Inequality(['E_A_72',   'ppErk'],  '<', ['E',        'ppErk'], 'ppErk_3')
        erk_ineq4 = Inequality(['E_A_72',   'ppErk'],  '<', ['D',        'ppErk'], 'ppErk_4')
        erk_ineq5 = Inequality(['E_A_72',   'ppErk'],  '<', ['T',        'ppErk'], 'ppErk_5')
        erk_ineq6 = Inequality(['A_72',     'ppErk'],  '<', ['E',        'ppErk'], 'ppErk_6')
        erk_ineq7 = Inequality(['A_72',     'ppErk'],  '<', ['D',        'ppErk'], 'ppErk_7')
        erk_ineq8 = Inequality(['A_72',     'ppErk'],  '<', ['T',        'ppErk'], 'ppErk_8')
        erk_ineq9 = Inequality(['E_M_72',   'ppErk'],  '>', ['E',        'ppErk'], 'ppErk_9')
        erk_ineq10 = Inequality(['M_72',     'ppErk'],  '<', ['E_M_72',   'ppErk'],'ppErk_10')

        erk = [
            erk_ineq1,
            erk_ineq2,
            erk_ineq3,
            erk_ineq4,
            erk_ineq5,
            erk_ineq6,
            erk_ineq7,
            erk_ineq8,
            erk_ineq9,
            erk_ineq10,
        ]
        s6k_ineq1 = Inequality(['E',        'pS6K'],  '<', ['T',        'pS6K'], 'pS6K_1')
        s6k_ineq2 = Inequality(['E',        'pS6K'],  '<', ['D',        'pS6K'], 'pS6K_2')
        s6k_ineq3 = Inequality(['E_A_72',   'pS6K'],  '<', ['T',        'pS6K'], 'pS6K_3')
        s6k_ineq4 = Inequality(['E_A_72',   'pS6K'],  '<', ['D',        'pS6K'], 'pS6K_4')
        s6k_ineq5 = Inequality(['E_A_72',   'pS6K'],  '<', ['A_72',     'pS6K'], 'pS6K_5')
        s6k_ineq6 = Inequality(['E_M_72',   'pS6K'],  '<', ['D',        'pS6K'], 'pS6K_6')
        s6k_ineq7 = Inequality(['E_M_72',   'pS6K'],  '<', ['T',        'pS6K'], 'pS6K_7')
        s6k_ineq8 = Inequality(['M_72',     'pS6K'],  '>', ['E_M_72',   'pS6K'], 'pS6K_8')

        s6k = [
            s6k_ineq1,
            s6k_ineq2,
            s6k_ineq3,
            s6k_ineq4,
            s6k_ineq5,
            s6k_ineq6,
            s6k_ineq7,
            s6k_ineq8,
        ]
        smad2_ineq1 = Inequality(['T',      'pSmad2'],  '>', ['D',           'pSmad2'], 'pSmad1_1')
        smad2_ineq2 = Inequality(['E',      'pSmad2'],  '>', ['T',           'pSmad2'], 'pSmad1_2')
        smad2_ineq3 = Inequality(['A_72',   'pSmad2'],  '>', ['A_1.25',      'pSmad2'], 'pSmad1_3')
        smad2_ineq4 = Inequality(['E_M_72', 'pSmad2'],  '<', ['E',           'pSmad2'], 'pSmad1_4')
        smad2_ineq5 = Inequality(['E_M_72', 'pSmad2'],  '<', ['E_M_1.25',    'pSmad2'], 'pSmad1_5')
        smad2_ineq6 = Inequality(['M_72',   'pSmad2'],  '<', ['E',           'pSmad2'], 'pSmad1_6')

        smad = [
            smad2_ineq1,
            smad2_ineq2,
            smad2_ineq3,
            smad2_ineq4,
            smad2_ineq5,
            smad2_ineq6,
        ]
        all_ineq = akt + erk + s6k + smad
        ineq = InequalityGroup(all_ineq)

        # SemiRSS()
        akt_mapping = Mapping('pAkt', 'pAkt')
        erk_mapping = Mapping('pErk', 'ppErk')
        s6k_mapping = Mapping('pS6K', 'pS6K')
        smad2_mapping = Mapping('pSmad2', 'pSmad2')

        mappings = [akt_mapping, erk_mapping, s6k_mapping, smad2_mapping]

        O = RandomSimulation(
            cross_talk_model_antstr(), exp_data=exp_data, inequalities= ineq, mappings=mappings,
            free_parameters=free_parameters, iterations=100
            )

        O.fit()


        # do_robustness(model_string)


    if OPTIMIZE:
        print(working_directory)
        azd_data = os.path.join(working_directory, r'data\HardCopy\AZD_calculations - v3.xlsx')
        mk_data = os.path.join(working_directory, r'data\HardCopy\MK2206_calculations - v3.xlsx')

        assert os.path.isfile(azd_data)
        assert os.path.isfile(mk_data)

        azd_data = pandas.read_excel(azd_data, sheet_name='AZDAverage', index_col='Condition')
        mk_data = pandas.read_excel(mk_data, sheet_name='MKAverage', index_col='Condition')

        mk_data = mk_data.iloc[3:]

        exp_data = pandas.concat([azd_data, mk_data])
        # print(exp_data)

        akt_mapping = Mapping('pAkt', 'pAkt')
        erk_mapping = Mapping('pErk', 'pErk')
        s6k_mapping = Mapping('pS6K', 'pS6K')
        smad2_mapping = Mapping('pSmad2', 'pSmad2')

        mappings = [akt_mapping, erk_mapping, s6k_mapping, smad2_mapping]

        free_parameters = OrderedDict({
            'kmTORC1Phos_ki': 0.001,
            'kPI3KPhosByTGFbR_kcat': 50.0,
            'kAktDephos_Vmax': 31.1252344504785,
            'kPI3KDephosByErk': 5.014,
            'kS6KPhosBymTORC1_kcat': 2.77975221288272,
            'kPI3KPhosByGF': 0.239474698704283,
            'kPI3KDephosByS6K': 25.0,
            'kErkPhos_kcat1': 85.0103161451182,
            'kmTORC1Dephos_Vmax': 1.0,
            'kS6KDephos_Vmax': 50.0,
            'kAktPhos_kcat': 2.9215,
            'kRafPhos_ki': 3.5,
            'kRafPhosByTGFbR_kcat': 265.0,
            'kRafPhosByPI3K_kcat': 50.0,
            'kMekPhos_kcat1': 149.5209856,
            'kMekPhos_ki1': 0.25,
            'kTGFbOn': 0.100647860357268,
            'kSmad2PhosByAkt_kcat': 1.0,
            'kSmad2Dephos_Vmax': 58.8712661228653,
            'kAktPhos_ki': 0.01,
            'kmTORC1Phos_kcat': 0.1,
        })
        free_parameters = free_parameters.keys()

        # o = RandomSimulation(cross_talk_model_antstr(), exp_data * 10, mappings, iterations=10,
        #                      free_parameters=free_parameters)
        # RSS = o._compute_obj_fun()
        # print(RSS.residuals)
        # print(RSS.rss)

        # print(exp_data*10)











