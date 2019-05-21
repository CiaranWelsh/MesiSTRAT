import PyBoolNet as pbn




string = '''
v1, !v1
v2, 1
v3, v2 & (!v1 | v3)
'''


primes = pbn.Repository.get_primes("remy_tumorigenesis")



print(primes)

print(pbn.FileExchange.primes2bnet(primes))

pbn.InteractionGraphs.create_empty_igraph(primes, "igraph.pdf")












