pub mod TrieStore {
    use trie_rs::{Trie, TrieBuilder};

    /// A trie that stores u64 values
    pub struct TrieStore {
        trie: Trie<u8>,
    }

    impl TrieStore {
        pub fn get(&self, key: &str) -> Option<u64> {
            self.trie
                .predictive_search(key.as_bytes())
                .first()
                .map(|bytes| {
                    let mut bytes = bytes.to_vec();
                    let value_bytes = bytes.split_off(bytes.len() - 8);
                    u64::from_be_bytes(value_bytes.try_into().unwrap())
                })
        }
    }

    pub struct TrieStoreBuilder {
        builder: TrieBuilder<u8>,
    }

    impl TrieStoreBuilder {
        pub fn new() -> Self {
            TrieStoreBuilder {
                builder: TrieBuilder::new(),
            }
        }

        pub fn push(&mut self, key: &str, value: u64) {
            let mut key = key.as_bytes().to_vec();
            key.extend_from_slice(&value.to_be_bytes());
            self.builder.push(key);
        }

        pub fn build(self) -> TrieStore {
            let trie = self.builder.build();
            TrieStore { trie }
        }
    }
}