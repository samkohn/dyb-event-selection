#include <vector>

TChain* get_chain(const char* chain_name, const char* path_format, std::vector<int> items)
{
    TChain* ch = new TChain(chain_name);
    char buffer[500];
    for(size_t i = 0; i < items.size(); i++)
    {
        sprintf(buffer, path_format, items[i]);
        ch->Add(buffer);
    }
    return ch;
}

