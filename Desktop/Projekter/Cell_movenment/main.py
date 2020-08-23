#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
plt.close('all')


def random_position():
    """  returns random point in plane (r1, r2, 0)  """
    pos = np.random.randn(3)
    pos[2] = 0
    return pos
    

# Class for cell
class Cell():
    """
    Contains information about the individual cell
    where it is, who it's connected to
    and what cell material it physically contains
    """
    def __init__(self, idx = 0, sugar=1, water=1, growth=0, connections=[], position=random_position(), force=np.zeros(3)):            
        # Outside cell
        self.position = position
        self.velocity = np.zeros(3)
        self.force = force
        self.connections = connections
        # Inside cell
        self.idx = idx
        self.sugar = sugar
        self.water = water
        self.growth = growth

"""
    def plot(self):
        #fig = plt.figure()
        #plt.clf()
        ax = plt.axes(projection='3d')
        
        for cell in self.cells:
            if cell.idx == 0:
                plot_style = 'or'
            elif cell.idx == 1:
                plot_style = 'ob'
            else:
                plot_style = 'ok'
                    
            x, y, z = cell.position[0], cell.position[1], cell.position[2]
            ax.plot3D([x], [y], [z], plot_style)
            ax.set_xlim(-2, 3)
            ax.set_ylim(-2, 3)
            ax.set_zlim(-2, 3)
            
        plt.pause(0.1)
"""

# Class for plant
class Plant():
    """
    Contains information about the plant as a whole
    what cells there are, how they interact
    as well as movement of the plant as a whole, 
    both regarding positions and chemicals of the cell. 
    """
    def __init__(self, cells = []):
        # Instantiate list of cells
        self.cells = self.create_seed() if len(cells) == 0 else cells
    
    
    def move_plant(self, dt = 0.1, drag = 0.1):
        """  
        Get forces applied to all cells in the plant 
        and update their velocity and position
        drag from 0-100% 
        """
        for cell in self.cells:
            cell.velocity *= 1 - drag # Apply drag
            cell.velocity += cell.force * dt # change speed [v(t+dt) = v(t) + f * dt because mass = 1 so f = a]
            cell.position += cell.velocity * dt # move cell [x(t+dt) = x(t) + v * dt]
        
        
    def get_forces(self, avg_dist=1, k=1):
            """  get forces and update forces of cells  """
            # Iterate all cells in plant
            for cell in self.cells:
                # Position of cell
                my_pos = cell.position
                # vector for summing forces
                force = np.zeros(3)
                
                # Iterate neighbours
                for neigh in cell.connections:
                    # Neighbours position
                    your_pos = self.cells[neigh].position
                    # Vector from you to me
                    me2you = my_pos - your_pos
                    # Vector from equilibrium to me
                    x = me2you * (1 - avg_dist / np.linalg.norm(me2you))
                    # Get hookes force
                    F = - k * x
                    # Sum of forces
                    force += F
                    
                # Update force vector of cell
                cell.force = force
            
            
    def create_seed(self, position_stem=np.array([0, 0, 0.]), position_root=np.array([0, 0, 2.])):
        """  Create initial stem and root cells that are connected to each other  """
        cells = []
        cell_stem = Cell(idx=0, position=position_stem, connections = [1])
        cell_root = Cell(idx=1, position=position_root, connections = [0])
        cells.append(cell_stem)
        cells.append(cell_root)            
        return cells
        
        
    def add_cell(self, connections=[], cell = None):
        """  Create a cell and add it to the plant  """
        # Create new cell
        idx = len(self.cells)
        position = np.random.randn(3) if cell == None else cell.position + 1
        connections = connections if len(connections) != 0 else [idx - 1]
        new_cell = Cell(idx=idx, connections=connections, position=position)
        # Connect new cell with connected cells in plant  """
        for old_cell in self.cells:
            if old_cell.idx in new_cell.connections: 
                old_cell.connections.append(new_cell.idx)
        # Add new cell to plant
        self.cells.append(new_cell)
        
        
    def position_matrix(self):
        """  Return all positions in an [N, 3] matrix  """
        positions = np.zeros([len(self.cells), 3])

        for i, cell in enumerate(self.cells):
            # Position of cell
            positions[i, :] = cell.position
        return positions
    
    
    def get_distances(self):
        """  Return distances between cells in an [N, N] matrix  """
        N = len(self.cells) # Number of cells
        distances = np.zeros([N, N]) # distances between cells
        positions = self.position_matrix() # positions of cells 
        
        # get distances between cells (exploit symmetry between upper and lower triangular form)
        for i, position in enumerate(positions[:-1, :]): # Iterate matrix except the last one
            directions = positions[i+1:, :] - position # direction from i to j > i
            distances[i, i+1:] = np.linalg.norm(directions, axis=1) # length of directions
        
        return distances + distances.T # Add lower triangle of matrix to upper      

        
    def get_connections(self):
        N = len(self.cells) # Number of cells
        connections = np.zeros([N, N]) # connections between cells
        
        for i, cell in enumerate(self.cells):
            connections[i, cell.connections] += 1
            
        return connections


    def connect_colliding_cells(self, collision_distance=1):
        # Get distances bigger than 0 (self interaction), and smaller than 1 (collision distance)
        too_close = (0 < self.get_distances()) * (self.get_distances() < collision_distance)
        not_connected = self.get_connections() == 0
        collisions = too_close * not_connected
        
        for x, y in np.array(np.nonzero(collisions)).T:
            print("Connected ",x, " and ", y)
            self.cells[x].connections.append(y)    

        
    def plot(self):
        # Create 3D axis
        ax = plt.axes(projection='3d')
        # 3d position where xyz[idx] = position_idx
        xyz = []
        for cell in self.cells:
            xyz.append(cell.position)
        xyz = np.array(xyz) # Turn into a [N, 3] numpy array
        
        # Get edges 
        edges = []
        for cell_idx, cell in enumerate(self.cells):
            for connected_cell_idx in cell.connections:
                if connected_cell_idx > cell_idx: # avoid double counting
                    # add two cell positions
                    edges.append(xyz[[cell_idx, connected_cell_idx], :])
        
        

            
        xs = xyz[:, 0]
        ys = xyz[:, 1]
        zs = xyz[:, 2]

        ax.plot3D(xs, ys, zs, 'og')
        
        for edge in edges:
            ax.plot3D(edge[:, 0], edge[:, 1], edge[:, 2], '-k')
        
        ax.set_xlim(-2, 2)
        ax.set_ylim(-2, 2)
        ax.set_zlim(-2, 2)
            
        plt.pause(0.1)        

        
class Environment():
    """
    Contains information about all factors in the environment: 
    plants, light, water, Mechanics, flow of nutrients etc. 
    Simulates the environment, and shows the simulation in progress. 
    """
    def __init__(self, plants=[]):
        # Make list of plants if none is given
        self.plants = self.create_plant() if len(plants) == 0 else plants
        
    def create_plant(self):
        # Make a plant and put it in a list
        return [Plant()]
    
    def move_plants(self, dt, plot_is_allowed = True):
        # 
        for plant in self.plants:
            plant.get_forces()
            plant.move_plant(dt=dt)
            if plot_is_allowed: plant.plot()
        
    def check_collision(self):
        for plant in self.plants: plant.connect_colliding_cells()
        
    def run(self, T = 5, dt = 0.2):
        # Run environment for T/dt steps of dt
        T_counter = 0.
        while T_counter <= T:
            self.move_plants(dt)
            self.check_collision()
            
            T_counter += dt

    
# -------------------------------------------------------------------------------------------------
# Initiate environment
my_env = Environment()
# Get first plant in environment
my_plant = my_env.plants[0]
# Get first cell in first plant
my_cell = my_plant.cells[0]

#%%
"""
Run simulation of environment
"""
# Time and timestep of simulation
T, dt = 0.5, 0.1

# Run environment, add plant, repeat 4 times
for _ in range(100):
    my_env.plants[0].add_cell()
    my_env.run(T=T, dt=dt)
    print("Added cell")


#%%
"""
Network - approach to plotting
Networkx Does not work in 3D :(
"""
"""
import matplotlib.pyplot as plt
import networkx as nx

plt.close('all')
plt.clf()
G = nx.Graph()



cells = my_env.plants[0].cells
for cell in cells:
    cell_idx = cell.idx
    cell_position = cell.position
    G.add_node(cell_idx, pos=list(cell_position[:2]))

for cell in cells:
    for connection in cell.connections:
        G.add_edge(cell.idx, cells[connection].idx)

pos = nx.get_node_attributes(G, 'pos')
nx.draw(G, pos)
"""
#%%

#%%
    
    
    

import networkx as nx

G=nx.Graph()

G.add_node(1,pos=(1,1))

G.add_node(2,pos=(2,2))

G.add_edge(1,2)

pos=nx.get_node_attributes(G,'pos')


nx.draw(G,pos)



nx.draw(G)
#nx.draw_networkx_nodes(G, pos, node_size=3000, nodelist=[0, 1, 2, 3, 4, 5], node_color='b')
#nx.draw_networkx_edges(G, pos, alpha=0.5, width=6)
plt.axis('off')
plt.show()


#%%
# explicitly set positions
pos = {0: (0, 0),
       1: (1, 0),
       2: (0, 1),
       3: (1, 1),
       4: (0.5, 2.0),
       5: (0, 2)}

#nx.draw_networkx_nodes(G, pos, node_size=2000, nodelist=[4])




 

#%%
"""  DISCARD PILE (that may still be useful)  """

"""    
    def plot_plant(self):
        x, y, z = [], [], []
        for cell in self.cells:
            x_new, y_new, z_new = cell.position[0], cell.position[1], cell.position[2]
            x.append(x_new)
            y.append(y_new)
            z.append(z_new)
        fig = plt.figure()
        ax = plt.axes(projection='3d')
        ax.plot3D(x, y, z, 'or')
"""