
class CombinatorialModel:


    def __init__(self,
                 problem_directory,
                 data_directory,
                 mutually_exclusive_reactions=[],
                 fit='1_1',
                 run_mode='slurm',
                 copy_number=33,
                 randomize_start_values=True,
                 overwrite_config_file=True,
                 method='particle_swarm',
                 population_size=75,
                 swarm_size=50,
                 iteration_limit=2000,
                 number_of_generations=500,
                 lower_bound=0.001,
                 upper_bound=10000,
                 use_best_parameters=False):

        self.mutually_exclusive_reactions = mutually_exclusive_reactions
        if self.mutually_exclusive_reactions is not None:
            if not isinstance(self.mutually_exclusive_reactions, list):
                raise TypeError('expecting list but got {}'.format(type(self.mutually_exclusive_reactions)))
            for i in self.mutually_exclusive_reactions:
                if not isinstance(i, tuple):
                    raise TypeError('expecting tuple but got {}'.format(type(self.mutually_exclusive_reactions)))

        self._topology = 0
        self.problem_directory = problem_directory
        if not os.path.isdir(self.problem_directory):
            os.makedirs(self.problem_directory)
        self.data_dir = data_directory
        if not os.path.isdir(self.data_dir):
            raise ValueError(f'{self.data_dir} is not a directory')

        self.fit = fit
        self.run_mode = run_mode
        self.copy_number = copy_number
        self.randomize_start_values = randomize_start_values
        self.overwrite_config_file = overwrite_config_file
        self.method = method
        self.population_size = population_size
        self.swarm_size = swarm_size
        self.iteration_limit = iteration_limit
        self.number_of_generations = number_of_generations
        self.lower_bound = lower_bound
        self.upper_bound = upper_bound
        self.use_best_parameters = use_best_parameters

        self.cps_file = os.path.join(self.topology_dir, 'Topology{}'.format(self.topology))

        # dict of reactions that vary with topologies and another dict with corresponding hypothesis names
        self.model_variant_reactions, self.topology_names = self._model_variant_reactions()


    def _model_variant_reactions(self):
        """
        Get all methods that begin with 'extension_hypothesis' and return their values in a dict[number] = reaction_string

        This assembles the reactions that are not in every model and will later be combinatorially combined with the
        core model.

        Returns:

        """
        hypothesis_reactions = []
        hypothesis_reaction_names = []
        for i in dir(self):
            if i.startswith('extension_hypothesis'):
                hypothesis_reactions.append(getattr(self, i)())
                hypothesis_reaction_names.append(i.replace('extension_hypothesis_', ''))

        dct = OrderedDict()
        names = OrderedDict()
        for i in range(len(hypothesis_reactions)):
            dct[i] = hypothesis_reactions[i]
            names[i] = hypothesis_reaction_names[i]
        return dct, names

    def __len__(self):
        """
        Subtract 1 for 0 indexed python
        :return:
        """
        return len(list(self._get_combinations()))




    def __iter__(self):

        while self.topology < len(self):
            yield self.topology
            self.topology += 1

    def __getitem__(self, item):
        if not isinstance(item, int):
            raise TypeError('"item" should be of type int. Got "{}" instead'.format(type(item)))

        self.topology = item
        # self.model_specific_reactions = self._assembel_model_reactions()[item]
        return self

    @property
    def topology(self):
        return self._topology

    @topology.setter
    def topology(self, new):
        assert isinstance(new, int)
        self._topology = new

