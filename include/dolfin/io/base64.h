#include <string>

auto base64_encode(unsigned char const* , unsigned int len) -> std::string;
auto base64_decode(std::string const& s) -> std::string;

