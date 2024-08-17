### Atom (Graph):

At the center of the atomphys package there is an atom. We decided to define atom as a graph data structure. This is because it is natural to think about an atom as a graph, where the states are the vertices, and the transitions are the edges of the graph. 

It allows us then to build upon great work that has been done on graphs in order to plot atoms, explore them, list all the properties, and their vertices, and edges. 


#### Loading of an atom:

But it would be extremaly tideous to build up such structure from scratch. Instead one can load an atom directly from NIST database, or from a custom json database. 

To load an atom from NIST, you can call from_nist() function. As a name you need to parse the element name of an atom you are interested. For instance if you are interested in Calcium atom/ion you need to call from_nist('Ca')/from_nist('Ca+').

> `atomphys.from_nist(name, _ureg)`
>
> Returns an atom with states and transitions found in NIST database. If NIST is not complete for your purposes, the missing transitions can be added manually. 

Alternatively one can use:

> `atomphys.from_json(name, _ureg)`
>
> Returns an atom with states and transitions found in the custom database. Database has to be a json file with a format given below

```
{"name": "Ca", 
"states": [{"configuration": "3p6.4s2", "term": "1S0", "energy": "0.0 Ry"}, {"configuration": "3p6.4s.4p", "term": "3P0", "energy": "0.1381289573728982 Ry"}, ... ], 
"transitions": [{"A": "2.74e+03 s^-1", "state_i": {"energy": "0.0 Ry", "term": "1S0"}, "state_f": {"energy": "0.138604292491823", "term": "3P1"}}, {"A": "4.89e+05 s^-1", "state_i": {"energy": "0.1381289573728982 Ry", "term": "3P0"}, "state_f": {"energy": "0.1853094352973106", "term": "3D1"}}...]
}
```
Above example 'Ca' database was taken from {Mills 2018}.

#### Basic Usage:
Upon importing your atom, it allows you to interact with it in a similar way you would interact with a databse / or a graph. It makes it easy to interact in an easy way with NIST database, or any other database that you would use. For the functionalities of the atom one should refer to the API reference, which should list all possible actions that can be performed. 

### States (Node) and Transitions (Edge):
Atom is formed by collection of states (nodes), connected via transitions (edges). Each State and Transition is a separate object with properties and possible functions associated with it. The representation of the atom as a graph could be used for the optimal repumping schemes, or optimal control schemes. 


### Electric Field

<span style="color:red">TO DO </span>