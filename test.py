import MDAnalysis
import numpy as np
import numpy.linalg as nl
import sys
import networkx as nx
import matplotlib.pyplot as plt

'''  ####################    GLOBAL VARIABLES   ################################## '''
traj_start = 1
traj_end   = 5000
step_size  = 10


class Cluster_Analysis:
    
    G = nx.Graph()
    '''  ###########################     DEFINE PARAMETERS HERE    ###############################  '''
    cutoff     = 5        # Distance in Angstrom adjacent monomers must be within to be considered "connected"
    surf       = 5        # Distance in Angstrom from surface monomers must be within to be considered "adsorbed"

    def __init__(self, topology, trajectory, res_id, num_res):
        self.res_id  = res_id
        self.num_res = num_res
        self.u = Universe(topology, trajectory)
        self.x = self.u.select_atoms('resname {}'.format(res_id))

    def plot_Aggregates(self, plot_list):
        timestep, bound_aggs, free_aggs, bound_mono, free_mono = zip(*plot_list)
        plt.plot(timestep, bound_aggs, linewidth = 0.5, label='Bound Aggregates')
#        plt.plot(timestep, bound_mono, linewidth = 0.5, label='Bound Monomers')
#        plt.plot(timestep, free_mono, linewidth = 0.5, label='Free Monomers')
        plt.plot(timestep, free_aggs, linewidth = 0.5, label='Free Aggregates')
        plt.legend()
        plt.show()

    def cluster_Classify(self, clust) -> str:
        # This seems rather arbitrary but it is retreiving the z component of two surface oxygens
        # from the top and bottom layers. These are then used to see if any clusters are adsorbed to the surface.
        # It assumes that the cartesian coordinates are non negative in the z component for both the top and
        # bottom layers and that the z component of the top > z component of the bot.
        # You may need to verify this in VMD and translate a trajectory if this is not the case.
        top   =  self.u.select_atoms('resname PYR and resnum 127 and name OB16').positions[0][2] - self.surf
        bot   =  self.u.select_atoms('resname PYR and resnum 188 and name OB11').positions[0][2] + self.surf
        
        for c in clust:
            temp_sel = self.u.select_atoms('resname {} and resnum {}'.format(self.res_id, c, updating='TRUE')).center_of_geometry()
            agg_size = len(clust)
            if agg_size == 1:
                if (temp_sel[2] < bot or temp_sel[2] > top):
                    return 'bound monomer', 1
                else:
                    return 'free monomer', 1
            else:
                if (temp_sel[2] < bot or temp_sel[2] > top):
                    return 'bound aggregate', agg_size
            return 'free aggregate', agg_size


    def get_Graph(self, edge_list):
        self.G.clear()
        self.G.add_edges_from(edge_list)


    def run(self, start, end, stepsize):
        plot_list = []
        for ts in self.u.trajectory[start:end:stepsize]:
            edge_list = []
            for residue in self.x.residues.resids:
                temp = self.u.select_atoms('resname {} and around {} resnum {}'.format(self.res_id, self.cutoff, residue), updating='TRUE')
                if (len(temp.residues.resids) == 0):
                    edge_list.append((residue, residue))
                else:
                    for node in temp.residues.resids:
                        edge_list.append((residue, node))

            # Builds a graph for each time step, sorts by cluster size
            self.get_Graph(edge_list)
            clusters = sorted(nx.connected_components(self.G), key=len, reverse=True)

            bound_aggregates, free_aggregates, bound_monomers, free_monomers = 0, 0, 0, 0
            for clust in clusters:
                output = self.cluster_Classify(clust)
                if (output[0] == 'bound aggregate'):
                    bound_aggregates = bound_aggregates + output[1]
                elif (output[0] == 'free aggregate'):
                    free_aggregates = free_aggregates + output[1]
                elif (output[0] == 'bound monomer'):
                    bound_monomers = bound_monomers + output[1]
                elif (output[0] == 'free monomer'):
                    free_monomers = free_monomers + output[1]
            plot_list.append((ts.time, bound_aggregates, free_aggregates, bound_monomers, free_monomers))
        self.plot_Aggregates(plot_list)
        


x = Cluster_Analysis(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
x.run(traj_start, traj_end, step_size)