import cobra

from SALib.sample import saltelli

import numpy as np
import scipy

import argparse
import pickle

try:
    from mpi4py import MPI
except:
    print("MPI not available: disabled")

if __name__ == '__main__':

    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()

    # Parse arguments.
    parser = argparse.ArgumentParser()
    parser.add_argument("model",
                        help="Path of constraint-based matlab model to analyze.")
    parser.add_argument("-N",
                        "--num_samples",
                        help="The number of samples to generate for each target reaction.",
                        type=int)
    args = parser.parse_args()
    
    # Load genome-wide constraint-based metabolic model.
    print("Rank {0} loading model {1}".format(rank, args.model))
    model = cobra.io.load_matlab_model(args.model)

    # Broadcast target reactions across tasks.
    if rank == 0:

        # Find inputs names.
        target = []
        for reaction in model.reactions:
            if 'EX_' in reaction.id and reaction.lower_bound < 0:
                target.append(reaction.id)

        print("This model contains:\n" +
              str(len(model.genes)) + " genes\n" +
              str(len(model.reactions)) + " reactions\n" +
              str(len(model.metabolites)) + " metabolites\n" +
              str(len(target)) + " target reactions\n")

    else:
        target = None

    target = comm.bcast(target, root=0)

    # Number of inputs.
    D = len(target)

    # Defining the model inputs.
    problem = {
        'num_vars' : D,
        'names' : target,
        'bounds' : [ [-10, 0] ] * D
    }

    # Number of samples to generate.
    N = args.num_samples
    if N is None:
        N = 2**4

    if rank == 0:

        # Return a NumPy matrix containing the model inputs using
        # Saltelli's sampling scheme.
        # Saltelli's scheme extends the Sobol sequence in a way to reduce
        # the error rates in the resulting sensitivity index calculations.
        # If calc_second_order is True
        # the resulting matrix has N * (2 * D + 2) rows
        # otherwise N * (D + 2) rows.

        X = saltelli.sample(problem,
                            N,
                            calc_second_order=False)

        print("The Saltelli sampler generated {0} samples.".format(X.shape[0]))

        # Model outputs.
        Y = np.empty(X.shape[0], dtype='float64')

        # Split the NumPy object by the number of available tasks.
        # Even if 'size' does not equally divide the axis.
        # For an array of length l that should be splitted into n sections,
        # it returns l % n sub-arrays of size l//n + 1 and the rest of size l//n.
        chunks = np.array_split(X, size)

        # Send one chunk per slave.
        count = np.empty(0)
        for i in range (1, len(chunks)):
            shape = chunks[i].shape
            comm.send(shape, dest=i, tag=i)
            comm.Send([chunks[i], MPI.DOUBLE], dest=i, tag=i+size)
            count = np.append(count, shape[0])

        # Master chunk.
        chunk = chunks[0]
        shape = chunks[0].shape

        # Add first dimension of master chunk.
        count = np.insert(count, 0, shape[0])

        # Cumulative sum of the elements of 'count'.
        # Insert a 0 value before index 0.
        # Return all elements except the last one.
        displ = np.insert(np.cumsum(count), 0, 0)[0:-1]

    else:
        X = None
        Y = None
        shape = None
        chunk = None
        count = None
        displ = None

    if rank != 0:
        shape = comm.recv(source=0, tag=rank)
        chunk = np.empty(shape, dtype='float64')
        print("Rank {0} received a chunk of shape {1} from rank 0.".format(rank, chunk.shape))
        comm.Recv([chunk, MPI.DOUBLE], source=0, tag=rank+size)

    # High Performance Computing (HPC) loop
    names = problem['names']
    get_reaction_by_id = model.reactions.get_by_id
    partial_Y = np.empty(chunk.shape[0], dtype='float64')

    print("Rank {0} is optimizing its chunk using flux balance analysis.".format(rank))
    for i, sample in enumerate(chunk):
        for name, value in zip(names, sample):
            get_reaction_by_id(name).lower_bound = value
            
            # print("sample of length = " + str(len(sample)) + " : " + str(sample) + "\n")
            # print("value : " + str(value) + "\n")
            # print("name : " + name + "\n")
            
        partial_Y[i] = model.slim_optimize()

    count = comm.bcast(count, root=0)
    displ = comm.bcast(displ, root=0)

    # Gather NumPy objects (following rank order).
    comm.Gatherv(partial_Y,
                 [Y, count, displ, MPI.DOUBLE],
                 root=0)

    # Global Synchronisation.
    comm.Barrier()

    # Write results.
    if rank == 0:

        with open('Problem.pkl', 'wb') as f:
            pickle.dump(problem, f, pickle.HIGHEST_PROTOCOL)

        with open('Y.pkl', 'wb') as f:
            pickle.dump(Y, f, pickle.HIGHEST_PROTOCOL)
