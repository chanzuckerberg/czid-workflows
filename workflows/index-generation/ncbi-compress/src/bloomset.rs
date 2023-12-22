use rayon::prelude::*;

pub struct BloomSet<const N: usize> {
    set: [u64; N],
    n_adds: u64,
}

impl<const N: usize> BloomSet<N> {
    // BloomSet operations are explicitly not parallelized, because the bloom set arrays are quite
    // small and the overhead of parallelization is not worth it.
    pub fn new() -> Self {
        BloomSet {
            set: [0; N],
            n_adds: 0,
        }
    }

    pub fn new_with_hashes(hashes: &[u64]) -> Self {
        let mut bloom_set = BloomSet::new();
        for hash in hashes {
            bloom_set.add(*hash);
        }
        bloom_set
    }

    pub fn add(&mut self, key: u64) {
        // Add to the count of items
        self.n_adds += 1;

        // Mod the hash by the total number of bits we have to get which bit we want to set
        let pos = key % (self.set.len() as u64 * 64);

        // This corresponds to which number contains the bit we want to set
        let which_number = pos / 64;
        // This corresponds to which bit in that number we want to set
        let which_bit = pos % 64;
        // Set the bit
        self.set[which_number as usize] |= 1 << which_bit;
    }

    pub fn union_mut(&mut self, other: &BloomSet<N>) {
        // OR each from other into self
        for (i, other) in other.set.iter().enumerate() {
            self.set[i] |= other;
        }
    }

    pub fn size(&self) -> u64 {
        self.set
            .iter()
            .fold(0, |acc, s| acc + (s.count_ones() as u64))
    }

    pub fn intersection_size(&self, other: &BloomSet<N>) -> u64 {
        self.set
            .iter()
            .zip(other.set)
            .map(|(a, b)| a & b)
            .fold(0, |acc, s| acc + (s.count_ones() as u64))
    }
}

pub trait BloomSetTreeable {
    fn hashes(&self) -> Vec<u64>;
    fn containment(&self, haystack: &Self) -> f64;
}

struct BloomSetTreeNode<T: BloomSetTreeable, const N: usize> {
    own: T,
    children_aggregate: BloomSet<N>,
}

pub struct BloomSetTree<T: BloomSetTreeable, const N: usize> {
    branching_factor: usize,
    nodes: Vec<BloomSetTreeNode<T, N>>,
}

impl<T: BloomSetTreeable + std::marker::Send + std::marker::Sync, const N: usize>
    BloomSetTree<T, N>
{
    pub fn new(branching_factor: usize) -> Self {
        BloomSetTree {
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

    pub fn insert(&mut self, item: T) {
        let item_bloom = BloomSet::new_with_hashes(&item.hashes());
        let node = BloomSetTreeNode {
            own: item,
            children_aggregate: item_bloom,
        };

        let mut current_idx = self.nodes.len();
        while let Some(parent_idx) = self.parent_idx(current_idx) {
            self.nodes[parent_idx]
                .children_aggregate
                .union_mut(&node.children_aggregate);
            current_idx = parent_idx;
        }
        self.nodes.push(node);
    }

    pub fn search(&self, item: &T, similarity_threshold: f64) -> Option<&T> {
        if self.nodes.is_empty() {
            return None;
        }

        let item_bloom = BloomSet::new_with_hashes(&item.hashes());
        let mut to_visit = vec![0];
        while !to_visit.is_empty() {
            let result = to_visit.par_iter().find_any(|node_idx| {
                let node = self.nodes.get(**node_idx).unwrap();
                item.containment(&node.own) >= similarity_threshold
            });

            if let Some(node_idx) = result {
                return Some(&self.nodes.get(*node_idx).unwrap().own);
            }

            to_visit = to_visit
                .par_iter()
                .flat_map(|node_idx| {
                    let node = self.nodes.get(*node_idx).unwrap();
                    if item_bloom.intersection_size(&node.children_aggregate) as f64
                        / item.hashes().len() as f64
                        >= similarity_threshold
                    {
                        self.child_idxes(*node_idx)
                    } else {
                        vec![]
                    }
                })
                .collect();
        }
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_intersection_size() {
        let mut a: BloomSet<64> = BloomSet::new();
        a.add(1);
        a.add(2);

        let mut b = BloomSet::new();
        b.add(2);
        b.add(3);

        assert_eq!(a.intersection_size(&b), 1);
    }
}
