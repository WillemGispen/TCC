#ifndef TCC_BCC15_H
#define TCC_BCC15_H

void Clusters_GetBCC_15();

int count_bonds_to_7A_ring(int first_7A_id, int new_particle_id);

void Cluster_Write_BCC15(int *first_7A_cluster, int new_particle_id);

#endif
