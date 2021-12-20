#ifndef SEQ_COMMON_H
#define SEQ_COMMON_H

//convert DNA (or RNA) ASCII characters to 2 bit numbers
//A->0, C->1, G-> 2, T->3, U->3
//a->0, c->1, g-> 2, t->3, u->3
unsigned char seq_nt4_table[256] = {
	0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

//From https://gist.github.com/badboy/6267743#64-bit-mix-functions
uint64_t hash64(uint64_t key, uint64_t mask)
{
  key = (~key + (key << 21)) & mask; // key = (key << 21) - key - 1;
  key = key ^ key >> 24;
  key = ((key + (key << 3)) + (key << 8)) & mask; // key * 265
  key = key ^ key >> 14;
  key = ((key + (key << 2)) + (key << 4)) & mask; // key * 21
  key = key ^ key >> 28;
  key = (key + (key << 31)) & mask;
  return key;
}

//From https://gist.github.com/badboy/6267743#32-bit-mix-functions
uint32_t hash32(uint32_t key, uint32_t mask)
{
  key = (~key + (key << 15)) & mask; // key = (key << 15) - key - 1;
  key = key ^ (key >> 12);
  key = (key + (key << 2)) & mask;
  key = key ^ (key >> 4);
  key = (key * 2057) & mask; // key = (key + (key << 3)) + (key << 11);
  key = key ^ (key >> 16);
  return key;
}

//Build inverse of unordered_map (assuming it is injective fn)
template<typename K, typename V>
void inverse_map (const std::unordered_map<K, V> &map, std::unordered_map<V, K> &inv)
{
  assert (inv.size() == 0);
  std::for_each(map.begin(), map.end(),
      [&inv] (const std::pair<K, V> &p) {
        assert (inv.find(p.second) == inv.end()); //key does not exist already
        inv.insert(std::make_pair(p.second, p.first));
      });
  assert (inv.size() == map.size());
}


#endif
