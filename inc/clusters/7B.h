#ifndef TCC_7B_H
#define TCC_7B_H

void Clusters_Get7B();

int count_bonds_to_6A_ring(int first_6A_id, int new_particle_id);
int count_double_bonds_to_6A_ring(int first_7A_id, int new_particle_id, int new_particle_id_);

void Cluster_Write_7B(int *first_6A_cluster, int new_particle_id);

#endif
