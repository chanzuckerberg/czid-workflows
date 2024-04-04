// this code is unused at this point, leaving in case we want to recover this later
// there was an issue with child aggregate hashes taking too much memory

pub mod minhashtree {
    use rayon::prelude::*;
    use sourmash::errors::SourmashError;
    use sourmash::sketch::minhash::KmerMinHash;

    use crate::ncbi_compress::ncbi_compress::containment;

    struct MinHashTreeNode {
        own: KmerMinHash,
        children_aggregate: KmerMinHash,
    }

    pub struct MinHashTree {
        branching_factor: usize,
        nodes: Vec<MinHashTreeNode>,
    }

    pub trait MinHashTreeFunctionality {
        fn new(branching_factor: usize) -> MinHashTree;
        fn parent_idx(&self, node: usize) -> Option<usize>;
        fn child_idxes(&self, node: usize) -> Vec<usize>;
        fn merge_to_parent(
            &mut self,
            parent_idx: usize,
            child_idx: usize,
        ) -> Result<(), SourmashError>;
        fn insert(&mut self, hash: KmerMinHash) -> Result<(), SourmashError>;
        fn contains(
            &self,
            hash: &KmerMinHash,
            similarity_threshold: f64,
        ) -> Result<bool, SourmashError>;
    }

    impl MinHashTreeFunctionality for MinHashTree {
        fn new(branching_factor: usize) -> Self {
            MinHashTree {
                branching_factor,
                nodes: Vec::new(),
            }
        }
        fn parent_idx(&self, node: usize) -> Option<usize> {
            if node == 0 {
                None
            } else {
                Some((node - 1) / self.branching_factor)
            }
        }
        fn child_idxes(&self, node: usize) -> Vec<usize> {
            let branching_factor = self.branching_factor;
            let first = branching_factor * node + 0 + 1;
            let last = branching_factor * node + branching_factor + 1;
            (first..last.min(self.nodes.len() - 1)).collect()
        }
        fn merge_to_parent(
            &mut self,
            parent_idx: usize,
            child_idx: usize,
        ) -> Result<(), SourmashError> {
            let (left, right) = self.nodes.split_at_mut(child_idx);

            left[parent_idx]
                .children_aggregate
                .merge(&right[0].children_aggregate)
        }
        fn insert(&mut self, hash: KmerMinHash) -> Result<(), SourmashError> {
            let node = MinHashTreeNode {
                own: hash.clone(),
                children_aggregate: hash.clone(),
            };
            let mut current_idx = self.nodes.len();
            self.nodes.push(node);

            while let Some(parent_idx) = self.parent_idx(current_idx) {
                self.merge_to_parent(parent_idx, current_idx)?;
                current_idx = parent_idx;
            }
            Ok(())
        }
        fn contains(
            &self,
            hash: &KmerMinHash,
            similarity_threshold: f64,
        ) -> Result<bool, SourmashError> {
            if self.nodes.is_empty() {
                return Ok(false);
            }

            let mut to_visit = vec![0];
            while !to_visit.is_empty() {
                let found = to_visit.par_iter().any(|node_idx| {
                    let node = self.nodes.get(*node_idx).unwrap();
                    containment(hash, &node.own).unwrap() >= similarity_threshold
                });

                if found {
                    return Ok(true);
                }

                to_visit = to_visit
                    .par_iter()
                    .flat_map(|node_idx| {
                        let node = self.nodes.get(*node_idx).unwrap();
                        if containment(hash, &node.children_aggregate).unwrap()
                            >= similarity_threshold
                        {
                            self.child_idxes(*node_idx)
                        } else {
                            vec![]
                        }
                    })
                    .collect();
            }
            Ok(false)
        }
    }
}
