import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import os, glob, pandas, numpy
from pycotools import model, viz, tasks
import tellurium as te
import seaborn
from collections import OrderedDict


seaborn.set_style('white')
seaborn.set_context('talk', font_scale=1)

"""
These are arguments for the conditions simulation functions. 
    :param GF: starting amount of GrowthFactors
    :param pretreatment: either 'AZD', or 'MK2206'. This is added as an event at the time specified by pretreatment time
    :param pretreatment_time: The time at which the variable specified by 'pretreatment' is added
    :param EV: Starting amount of EVerolimus
    :param serum_starve_event: Boolean, whether to remove serum, aka GrowthFactors by event
    :param TGFb_event: Boolean. Whether to add 1 unit TGF at t=-45min (71.75h)
"""

AZD_CONDITIONS = OrderedDict()
AZD_CONDITIONS['D']             =   [1,  None,     None,       0,  True,   False]
AZD_CONDITIONS['T']             =   [1,  None,     None,       0,  True,   True]
AZD_CONDITIONS['T_E']           =   [1,  None,     None,       1,  True,   True]
AZD_CONDITIONS['T_A_E_0']       =   [1,  "AZD",    0.0,        1,  True,   True]
AZD_CONDITIONS['T_A_E_24']      =   [1,  "AZD",    24,         1,  True,   True]
AZD_CONDITIONS['T_A_E_48']      =   [1,  "AZD",    48.0,       1,  True,   True]
AZD_CONDITIONS['T_A_E_70.75']   =   [1,  "AZD",    70.75,      1,  True,   True]
AZD_CONDITIONS['T_A_0']         =   [1,  "AZD",    0.0,        0,  True,   True]
AZD_CONDITIONS['T_A_24']        =   [1,  "AZD",    24,         0,  True,   True]
AZD_CONDITIONS['T_A_48']        =   [1,  "AZD",    48.0,       0,  True,   True]
AZD_CONDITIONS['T_A_70.75']     =   [1,  "AZD",    70.75,      0,  True,   True]

MK_CONDITIONS = OrderedDict()
MK_CONDITIONS['D']            =     [1, 0, None,      None,   0, True, False]
MK_CONDITIONS['T']            =     [1, 0, None,      None,   0, True, True]
MK_CONDITIONS['T_E']          =     [1, 0, None,      None,   1, True, True]
MK_CONDITIONS['T_M_E_0']      =     [1, 0, "MK2206",  0.0,    1, True, True]
MK_CONDITIONS['T_M_E_24']     =     [1, 0, "MK2206",  24.0,   1, True, True]
MK_CONDITIONS['T_M_E_48']     =     [1, 0, "MK2206",  48.0,   1, True, True]
MK_CONDITIONS['T_M_E_70.75']  =     [1, 0, "MK2206",  70.75,  1, True, True]
MK_CONDITIONS['T_M_0']        =     [1, 0, "MK2206",  0.0,    0, True, True]
MK_CONDITIONS['T_M_24']       =     [1, 0, "MK2206",  24.0,   0, True, True]
MK_CONDITIONS['T_M_48']       =     [1, 0, "MK2206",  48.0,   0, True, True]
MK_CONDITIONS['T_M_70.75']    =     [1, 0, "MK2206",  70.75,  0, True, True]


MODEL_SPECIES = ['TGFb', 'TGFbR', 'TGFbR_a', 'Smad2',
                 'pSmad2', 'Mek', 'pMek',
                 'Erk', 'pErk', 'GFR', 'pGFR', 'PI3K',
                 'pPI3K', 'Akt', 'pAkt', 'mTORC1',
                 'pmTORC1', 'S6K', 'pS6K']

MODEL_SPECIES += ['Everolimus', 'MK2206', 'AZD', 'GrowthFactors']

ALL_CONDITIONS = ['D', 'T', 'T_E',
                  'T_M_E_0', 'T_M_E_24',
                  'T_M_E_48', 'T_M_E_70.75',
                  'T_M_0', 'T_M_24', 'T_M_48',
                  'T_M_70.75', 'T_A_E_0', 'T_A_E_24',
                  'T_A_E_48', 'T_A_E_70.75', 'T_A_0',
                  'T_A_24', 'T_A_48', 'T_A_70.75']

def cross_talk_model_antstr():
    """
    Antimony string model containing the most currect version of
    the CrossTalk model from Katrine and Patricia's project.
    :return:
    """
    return """
    model TGFbModule()
        compartment Cell = 1.0
        
        var TGFbR       in Cell  
        var TGFbR_a     in Cell  
        var Smad2       in Cell  
        var pSmad2      in Cell  
        var Mek         in Cell
        var pMek        in Cell  
        var Erk         in Cell
        var pErk        in Cell  
        var GFR         in Cell  
        var pGFR        in Cell  
        var PI3K        in Cell  
        var pPI3K       in Cell  
        var Akt         in Cell
        var pAkt        in Cell  
        var mTORC1      in Cell  
        var pmTORC1     in Cell  
        var S6K         in Cell
        var pS6K        in Cell  
        
        
        var TGFb             in Cell
        var AZD              in Cell
        var GrowthFactors    in Cell
        var MK2206           in Cell
        var Everolimus       in Cell
    
        
        //TGFb module
        R1_1: TGFbR + TGFb    => TGFbR_a          ; Cell * kTGFbOn         *TGFb                ;
        R1_2: TGFbR           => TGFbR_a          ; Cell * kTGFbOnBasal    *TGFb                     ;
        R2:   TGFbR_a         => TGFbR + TGFb     ; Cell * kTGFbOff        *TGFbR_a             ;
        R3:   TGFb            =>                  ; Cell * kTGFRemoval     *TGFb                ;
        R4:   Smad2           => pSmad2           ; Cell * kSmad2Phos      *Smad2   *TGFbR_a    ;
        R5:   pSmad2          => Smad2            ; Cell * kSmad2Dephos    *pSmad2              ;
        //Rx1: Smad2          => pSmad2           ; Cell * kSmad2PhosByAkt * Smad2  *pAkt
        //kSmad2PhosByAkt = 10
        
        //MEK module                   //
        R9_1:    Mek      => pMek     ;   Cell * kMekPhosByPI3K      *Mek    *pPI3K   ;
        R9_2:    Mek      => pMek     ;   Cell * kMekPhosByTGFbR_a   *Mek    *TGFbR_a ;
        R9_3:    Mek      => pMek     ;   Cell * kMekPhosByGFR       *Mek    *pGFR    ;
        R10_1:   pMek     => Mek      ;   Cell * kMekDephosByAkt     *pMek   *pAkt    ;
        R10_2:   pMek     => Mek      ;   Cell * kMekDephosByAZD     *pMek   *AZD     ; 
        R11 :    Erk      => pErk     ;   Cell * kErkPhosByMek       *Erk    *pMek    ;
        R12:     pErk     => Erk      ;   Cell * kErkDephos          *pErk
        
        
        //PI3K Module
        R13:        GFR     => pGFR     ;   Cell *(kGFRPhos_kcat       * GrowthFactors * GFR /(kGFRPhos_km^kGFRPhos_h + GFR^kGFRPhos_h))
        R14:        pGFR    => GFR      ;   Cell *  kGFRDephos          *pGFR                    ;
        R15_1:      PI3K    => pPI3K    ;   Cell *  kPI3KPhosByGFR      *PI3K       *pGFR        ;
        R15_2:      PI3K    => pPI3K    ;   Cell *  kPI3KPhosByMek      *PI3K       *pMek        ;
        R15_1:      pPI3K   => PI3K     ;   Cell *  kPI3KDephosByS6K    *pPI3K      *pS6K        ;
        R16_2:      pPI3K   => PI3K     ;   Cell *  kPI3KDephosByErk    *pPI3K      *pErk        ;
        R17:        Akt     => pAkt     ;   Cell *  kAktPhos            *Akt        *pPI3K       ;
        R18_1:      pAkt    => Akt      ;   Cell *  kAktDephos          *pAkt                    ;
        R18_2:      pAkt    => Akt      ;   Cell *  kAktDephosByMK      *pAkt       *MK2206      ;
        R19_1:      mTORC1  => pmTORC1  ;   Cell *  kmTORC1PhosByAkt    *mTORC1     *pAkt        ;
        R19_2:      mTORC1  => pmTORC1  ;   Cell *  kmTORC1PhosByErk    *mTORC1     *pErk        ;
        R20_1:      pmTORC1 => mTORC1   ;   Cell *  kmTORC1DephosByEv   *pmTORC1    *Everolimus  ;
        R20_2:      pmTORC1 => mTORC1   ;   Cell *  kmTORC1Dephos       *pmTORC1                 ;
        R21_1:      S6K     => pS6K     ;   Cell *  kS6KPhosBymTORC1    *S6K        *pmTORC1     ;
        R21_2:      S6K     => pS6K     ;   Cell *  kS6KPhosByErk       *S6K        *pErk        ;
        R22:        pS6K    => S6K      ;   Cell *  kS6KDephos          *pS6K                    ;
        
        kTGFbOn             = 0.04
        kTGFbOnBasal        = 0.00004
        kTGFbOff            = 0.02
        kTGFRemoval         = 0.1
        kSmad2Phos          = 0.1
        kSmad2Dephos        = 0.01
          
        kMekPhosByPI3K      = 0.01    
        kMekPhosByTGFbR_a   = 0.01        
        kMekPhosByGFR       = 0.01    
        kMekDephosByAkt     = 0.01    
        kMekDephosByAZD     = 20  
        kErkPhosByMek       = 0.1
        kErkDephos          = 0.01
        
        kGFRPhos_kcat      = 10
        kGFRPhos_km        = 25
        kGFRPhos_h         = 2
        kGFRDephos         = 1
        kPI3KPhosByGFR     = 0.01    
        kPI3KPhosByMek      = 0.01    
        kPI3KDephosByS6K    = 0.1        
        kPI3KDephosByErk    = 0.001      
        kAktPhos            = 0.1
        kAktDephos          = 0.01
        kAktDephosByMK      = 0.5
        kmTORC1PhosByAkt    = 0.001
                 
        kmTORC1PhosByErk    = 0.00        // set to 0 on 05-12-2018 at 12:20
        kmTORC1Dephos       = 0.001   
        kmTORC1DephosByEv   = 1 
        kS6KPhosBymTORC1    = 0.1        
        kS6KPhosByErk       = 0.00        // set to 0 on 05-12-2018 at 12:20
        kS6KDephos          = 0.01
        
        //collect the feedback parameters into one list. Turn them off to see what they do
    
        TGFb                    = 0.0004  // set to 4e-4 on 05-12-2018 at 12:36 to include basal flux through TGF module
        AZD                     = 0
        GrowthFactors           = 1
        MK2206                  = 0
        Everolimus              = 0
        
        TGFbR                   = 100
        TGFbR_a                 = 0
        Smad2                   = 100
        pSmad2                  = 0
        Mek                     = 100
        pMek                    = 0
        Erk                     = 100
        pErk                    = 0
        GFR                     = 100
        pGFR                    = 0
        PI3K                    = 100
        pPI3K                   = 0
        Akt                     = 100
        pAkt                    = 0
        mTORC1                  = 100
        pmTORC1                 = 0
        S6K                     = 100
        pS6K                    = 0
    
      unit volume = 1 litre;
      unit time_unit = 3600 second;
      unit substance = 1e-9 mole;
        
    end;
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


def dose_response_for_growth_factors():
    """
    Use the dose repsonse function to see the effect of changing the
    initial concentration of GrowthFactors. Plot some output
    :return:
    """
    model_species = ['TGFbR', 'TGFbR_a',
                     'Smad2', 'pSmad2',
                     'Mek', 'pMek', 'Erk',
                     'pErk', 'GFR', 'pGFR',
                     'PI3K', 'pPI3K', 'Akt',
                     'pAkt', 'mTORC1', 'pmTORC1',
                     'S6K', 'pS6K', 'TGFb']
    mod = te.loada(cross_talk_model_antstr())

    l = [
        ['PI3K', 'pPI3K'],
        ['Mek', 'pMek'],
        ['Erk', 'pErk'],
        ['Akt', 'pAkt'],
        ['S6K', 'pS6K'],
        ['GFR', 'pGFR']
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


def add_serum_starve_event(model_string):
    """
    Set growth factors to 0 at t=70.25 minutes
    :param model_string:
    :return:
    """
    event_str = "RemoveGrowthFactors: at (time>70.25): GrowthFactors=0;\n"
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


def simulate_condition(model_string, GF=0, pretreatment=None,
                       pretreatment_time=None, EV=0,
                       serum_starve_event=False, TGFb_event=False, open_with_copasi=False):
    """
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
    if serum_starve_event:
        model_string = add_serum_starve_event(model_string)

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

    if open_with_copasi:
        ## has to be the modified model NOT model string
        copasi_mod = load_model_with_pyco(mod.getCurrentAntimony(), copasi_filename)
        tasks.TimeCourse(copasi_mod, start=0, end=72, step_size=1, intervals=73, run=False)

        copasi_mod.open()

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



def simulate_conditions_and_plot_as_bargraph(y, type='AZD'):
    """
    Takes output from simulation and plot
    :param df:
    :return:
    """
    seaborn.set_context('talk', font_scale=2)
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
        df = simulate_condition(cross_talk_model_antstr(), *conditions[k])
        df = df[df['time'] == 72.0]
        if df.empty:
            raise ValueError("Condition '{}' with values '{}' produces an "
                             "empty data frame".format(k, v))

        dct[k] = df

    df = pandas.concat(dct)

    df.index = df.index.droplevel(1)
    df = df.drop('time', axis=1)
    df = df.reset_index()
    df = df.rename(columns={'index': 'Condition'})

    if y not in df.columns:
        raise ValueError('The variable "{}" does not appear in your model'.format(y))

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
    print('Figure saved to "{}"'.format(fname))



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
        df = simulate_condition(cross_talk_model_antstr(), *conditions[k])
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


def make_condition(condition):
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

        model_string = cross_talk_model_antstr()

        try:
            GF, TGF, pretreatment, pretreatment_time, EV, serum_starve_event, TGFb_event = AZD_CONDITIONS[condition]
        except KeyError:
            GF, TGF, pretreatment, pretreatment_time, EV, serum_starve_event, TGFb_event = MK_CONDITIONS[condition]

        if serum_starve_event:
            model_string = add_serum_starve_event(model_string)

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
        mod.model['init([TGFb])'] = TGF
        mod.Everolimus = EV
        mod.GrowthFactors = GF

        if pretreatment_time == 0:
            if pretreatment == 'AZD':
                mod.AZD = 1
            elif pretreatment == 'MK2206':
                mod.MK2206 = 1

        return mod.getCurrentAntimony()

def open_condition_with_copasi(condition):
    """
    Create condition with model string and open with copasi,
    saving to folder while your at it
    :param condition:
    :return:
    """
    mod_string = make_condition(condition)
    copasi_dir = os.path.join(working_directory, 'copasi_models')
    if not os.path.isdir(copasi_dir):
        os.makedirs(copasi_dir)
    fname = os.path.join(copasi_dir, "{}.cps".format(condition))
    copasi_mod = load_model_with_pyco(mod_string, fname)
    tasks.TimeCourse(copasi_mod, start=0, end=72, step_size=1, intervals=73, run=False)
    copasi_mod.open()


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
            raise ValueError

    if kwargs.get('figsize') is not None:
        figsize = kwargs['figsize']
        del kwargs['figsize']
    else:
        figsize = (2*len(cond), 3*len(cond))

    cols = seaborn.color_palette("hls", len(vars))*len(cond)
    cols = iter(cols)

    fig = plt.figure(figsize=figsize)
    for i, c in enumerate(cond):
        if c not in MK_CONDITIONS.keys() + AZD_CONDITIONS.keys():
            raise ValueError('condition "{}" not in "{}"'.format(c, MK_CONDITIONS.keys() + AZD_CONDITIONS.keys()))

        model_str = make_condition(c)
        mod = te.loada(model_str)
        res = mod.simulate(0, 75, 750, ['time'] + vars)
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




if __name__ == '__main__':

    working_directory = r'D:\MesiSTRAT\CrossTalkModel'
    copasi_filename = os.path.join(working_directory, 'KatrinesTopology.cps')
    graphs_directory = os.path.join(working_directory, 'SimulationGraphs')

    """
    Set flags to determine which part of the script will run
    """

    DOSE_RESPONSE_GROWTH_FACTOR     = False
    DOSE_RESPONSE_TGFB              = False
    GET_ODES_WITH_ANTIMONY          = False
    GET_MODEL_AS_SBML               = False

    SIMULATE_TIME_SERIES            = False
    SIMULATE_BAR_GRAPHS             = False

    OPEN_CONDITION_WITH_COPASI      = True

    SIMULATE_INPUTS                 = False



    if not os.path.isdir(graphs_directory):
        os.makedirs(graphs_directory)

    if DOSE_RESPONSE_GROWTH_FACTOR:
        # i = 'GFR'
        # j = 'pGFR'
        # pmid_for_IGF_data = 'PMID: 26217307'
        # fig = dose_response(cross_talk_model_antstr(), 'GrowthFactors', 0.1, 1000, 50, [i, j])
        # plt.title('Dose Response of {} and \n{} to GrowthFactors'.format(i, j))
        # fname = os.path.join(graphs_directory, 'GrowthFactorsDoseResponse{}.png'.format(i))
        # fig.savefig(fname, dpi=150, bbox_inches='tight')
        # plt.show()
        dose_response_for_growth_factors()
        mod = load_model_with_pyco(cross_talk_model_antstr(), copasi_filename)
        mod.open()

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


    if SIMULATE_BAR_GRAPHS:
        simulate_all_conditions_and_plot_as_bargraphs()

    if OPEN_CONDITION_WITH_COPASI:
        open_condition_with_copasi('T_E')

    ## use model checking to evaluate the truth of a condition for each parameter value

    if SIMULATE_TIME_SERIES:
        simulate_model_component_timecourse(['pAkt'], AZD_CONDITIONS.keys(), filename='AZD_pAkt')
        simulate_model_component_timecourse(['pAkt'], MK_CONDITIONS.keys(), filename='MK_pAkt')

    if SIMULATE_BAR_GRAPHS:
        # simulate_conditions_and_plot_as_bargraph('pAkt', 'AZD')
        # simulate_conditions_and_plot_as_bargraph('pAkt', 'MK2206')
        simulate_all_conditions_and_plot_as_bargraphs()


    if SIMULATE_INPUTS:
        simulate_inputs_only('AZD')
        simulate_inputs_only('MK2206')




























