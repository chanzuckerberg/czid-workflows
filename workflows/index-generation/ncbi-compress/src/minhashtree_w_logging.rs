pub mod minhashtree_w_logging {
    use rayon::prelude::*;
    use sourmash::encodings::HashFunctions;
    use sourmash::errors::SourmashError;
    use sourmash::signature::SigsTrait;
    use sourmash::sketch::minhash::KmerMinHash;

    use crate::ncbi_compress::ncbi_compress::containment;

    struct LoggingMinHashTreeNode {
        // save accession to node if we have matches, useful for logging sequence retention
        own: KmerMinHash,
        children_aggregate: KmerMinHash,
        accession: String,
    }

    pub struct MinHashTreeWithLogging {
        branching_factor: usize,
        nodes: Vec<LoggingMinHashTreeNode>, 
    }

    pub trait MinHashTreeWithLoggingFunctionality {
        fn new(branching_factor: usize) -> MinHashTreeWithLogging;
        fn parent_idx(&self, node: usize) -> Option<usize>;
        fn child_idxes(&self, node: usize) -> Vec<usize>;
        fn merge_to_parent(
            &mut self,
            parent_idx: usize,
            child_idx: usize,
        ) -> Result<(), SourmashError>;
        fn insert(&mut self, hash: KmerMinHash, accession: &str) -> Result<(), SourmashError>;
        fn contains(
            &self,
            hash: &KmerMinHash,
            similarity_threshold: f64,
        ) -> Result<Option<Vec<(String, String)>>, SourmashError>;
    }

    impl MinHashTreeWithLoggingFunctionality for MinHashTreeWithLogging {
        fn new(branching_factor: usize) -> Self {
            MinHashTreeWithLogging {
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
        fn insert(&mut self, hash: KmerMinHash, accession: &str) -> Result<(), SourmashError> {
            let node = LoggingMinHashTreeNode {
                own: hash.clone(),
                children_aggregate: hash.clone(),
                accession: accession.to_string(),
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
        ) -> Result<Option<Vec<(String, String)>>, SourmashError> {
            if self.nodes.is_empty() {
                return Ok(None);
            }

            let mut to_visit = vec![0]; // was initially zero
            while !to_visit.is_empty() {
                // Check if any of the nodes in the to_visit list are similar enough
                let found_accession = to_visit.par_iter().filter_map(|node_idx| {
                    let node = self.nodes.get(*node_idx).unwrap();
                    let containment_value = containment(hash, &node.own).unwrap();
                    if containment_value >= similarity_threshold {
                        Some((node.accession.clone(), containment_value.to_string()))
                    } else {
                        None
                    }
                }).collect::<Vec<_>>();
                // If we found a similar node, we can stop searching
                if !found_accession.is_empty() {
                    return Ok(Some(found_accession));
                }
                // Otherwise, we need to search the children of the nodes in the to_visit list
                to_visit = to_visit
                    .par_iter()
                    .flat_map(|node_idx| {
                        let node = self.nodes.get(*node_idx).unwrap();
                        // If the children are similar enough, we need to search them
                        if containment(hash, &node.children_aggregate).unwrap() >= similarity_threshold
                        {
                            self.child_idxes(*node_idx)
                        } else {
                            vec![]
                        }
                    })
                    .collect();
            }
            // we didn't find any similar nodes in the tree
            Ok(None)
        }
    }
}