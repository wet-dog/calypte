#include <array>
#include <iostream>
#include <limits>

namespace {
enum class Material {
    KERATIN_CORTEX,
    KERATIN_SEPARATION,
    MELANIN_MEMBRANE,
    MELANIN_INTERNAL,
    MELANIN_INTERNAL_TOP,
    AIR
};

constexpr size_t num_layers = 51;

using layers_t    = std::array<Material, num_layers>;
using distances_t = std::array<uint32_t, num_layers>;
} // namespace

constexpr std::array<Material, num_layers> init_layers() {
    layers_t result = {};

    // Initialise layers
    size_t i    = 0;
    result[i++] = Material::AIR;
    result[i++] = Material::KERATIN_CORTEX;
    result[i++] = Material::MELANIN_MEMBRANE;
    result[i++] = Material::MELANIN_INTERNAL_TOP;
    result[i++] = Material::MELANIN_MEMBRANE;

    size_t n = 0;
    // 12 melanosome layers minus 1 because of the unique top melanosome layer
    size_t num_melanosome_layers = 11;
    while (n < num_melanosome_layers) {
        result[i++] = Material::KERATIN_SEPARATION;
        result[i++] = Material::MELANIN_MEMBRANE;
        result[i++] = Material::MELANIN_INTERNAL;
        result[i++] = Material::MELANIN_MEMBRANE;
        ++n;
    }

    result[i++] = Material::KERATIN_CORTEX;
    result[i++] = Material::AIR;

    return result;
}

constexpr distances_t init_d_list(layers_t layers) {
    distances_t d_list = {};

    // Thicknesses in nm
    constexpr uint32_t keratin_cortex_thickness      = 5;
    constexpr uint32_t keratin_separation_thickness  = 50;
    constexpr uint32_t melanosome_membrane_thickness = 30;
    // Note: this can be changed to adjust periodicity
    constexpr uint32_t melanosome_internal_thickness     = 100;
    constexpr uint32_t top_melanosome_internal_thickness = 50;

    size_t i = 0;
    for (const auto &material : layers) {
        switch (material) {
            case Material::KERATIN_CORTEX:
                d_list[i] = keratin_cortex_thickness;
                break;
            case Material::KERATIN_SEPARATION:
                d_list[i] = keratin_separation_thickness;
                break;
            case Material::MELANIN_MEMBRANE:
                d_list[i] = melanosome_membrane_thickness;
                break;
            case Material::MELANIN_INTERNAL:
                d_list[i] = melanosome_internal_thickness;
                break;
            case Material::MELANIN_INTERNAL_TOP:
                d_list[i] = top_melanosome_internal_thickness;
                break;
            case Material::AIR:
                // "Infinity"
                d_list[i] = std::numeric_limits<uint32_t>::max();
                break;
        }
        ++i;
    }

    return d_list;
}

template <std::size_t S>
constexpr std::array<uint32_t, S> init_n_list(Material material,
                                              layers_t layers) {
    std::array<uint32_t, S> n_list = {};

    size_t i = 0;
    size_t j = 0;
    for (const auto &m : layers) {
        if (m == material) {
            n_list[j++] = i;
        }
        ++i;
    }

    return n_list;
}

namespace {
constexpr layers_t layers = init_layers();

auto count = [](Material material) constexpr {
    size_t i = 0;
    for (const auto &m : layers) {
        if (m == material)
            ++i;
    }
    return i;
};

constexpr auto keratin_cortex_indices =
    init_n_list<count(Material::KERATIN_CORTEX)>(Material::KERATIN_CORTEX,
                                                 layers);

constexpr auto keratin_separation_indices =
    init_n_list<count(Material::KERATIN_SEPARATION)>(
        Material::KERATIN_SEPARATION, layers);

constexpr auto melanin_membrane_indices =
    init_n_list<count(Material::MELANIN_MEMBRANE)>(Material::MELANIN_MEMBRANE,
                                                   layers);

constexpr auto melanin_internal_indices =
    init_n_list<count(Material::MELANIN_INTERNAL)>(Material::MELANIN_INTERNAL,
                                                   layers);

constexpr auto melanin_internal_top_indices =
    init_n_list<count(Material::MELANIN_INTERNAL_TOP)>(
        Material::MELANIN_INTERNAL_TOP, layers);

constexpr auto air_indices =
    init_n_list<count(Material::AIR)>(Material::AIR, layers);

constexpr distances_t d_list = init_d_list(layers);
} // namespace

int main() {
    // for (size_t i = 0; i < num_layers; ++i) {
    //     std::cout << static_cast<size_t>(layers[i]) << '\n';
    // }

    std::cout << '\n';

    for (size_t i = 0; i < count(Material::KERATIN_CORTEX); ++i) {
        std::cout << keratin_cortex_indices[i] << ' ';
    }

    std::cout << '\n';

    for (size_t i = 0; i < count(Material::KERATIN_SEPARATION); ++i) {
        std::cout << keratin_separation_indices[i] << ' ';
    }

    std::cout << '\n';

    for (size_t i = 0; i < count(Material::MELANIN_MEMBRANE); ++i) {
        std::cout << melanin_membrane_indices[i] << ' ';
    }

    std::cout << '\n';

    for (size_t i = 0; i < num_layers; ++i) {
        std::cout << d_list[i] << '\n';
    }
}