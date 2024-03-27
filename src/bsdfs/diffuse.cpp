#include <mitsuba/core/properties.h>
#include <mitsuba/core/spectrum.h>
#include <mitsuba/core/warp.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/texture.h>

#include <Eigen/Core>

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

constexpr size_t k_cortex_count = count(Material::KERATIN_CORTEX);
constexpr auto k_cortex_indices_ =
    init_n_list<k_cortex_count>(Material::KERATIN_CORTEX, layers);

constexpr size_t k_separation_count = count(Material::KERATIN_SEPARATION);
constexpr auto k_separation_indices_ =
    init_n_list<k_separation_count>(Material::KERATIN_SEPARATION, layers);

constexpr size_t m_membrane_count = count(Material::MELANIN_MEMBRANE);
constexpr auto m_membrane_indices_ =
    init_n_list<m_membrane_count>(Material::MELANIN_MEMBRANE, layers);

constexpr size_t m_internal_count = count(Material::MELANIN_INTERNAL);
constexpr auto m_internal_indices_ =
    init_n_list<m_internal_count>(Material::MELANIN_INTERNAL, layers);

constexpr size_t m_internal_top_count = count(Material::MELANIN_INTERNAL_TOP);
constexpr auto m_internal_top_indices_ =
    init_n_list<m_internal_top_count>(Material::MELANIN_INTERNAL_TOP, layers);

constexpr size_t air_count  = count(Material::AIR);
constexpr auto air_indices_ = init_n_list<air_count>(Material::AIR, layers);

constexpr distances_t d_list = init_d_list(layers);
} // namespace

NAMESPACE_BEGIN(mitsuba)

/**!

.. _bsdf-diffuse:

Smooth diffuse material (:monosp:`diffuse`)
-------------------------------------------

.. pluginparameters::

 * - reflectance
   - |spectrum| or |texture|
   - Specifies the diffuse albedo of the material (Default: 0.5)

The smooth diffuse material (also referred to as *Lambertian*)
represents an ideally diffuse material with a user-specified amount of
reflectance. Any received illumination is scattered so that the surface
looks the same independently of the direction of observation.

.. subfigstart::
.. subfigure:: ../../resources/data/docs/images/render/bsdf_diffuse_plain.jpg
   :caption: Homogeneous reflectance
.. subfigure:: ../../resources/data/docs/images/render/bsdf_diffuse_textured.jpg
   :caption: Textured reflectance
.. subfigend::
   :label: fig-diffuse

Apart from a homogeneous reflectance value, the plugin can also accept
a nested or referenced texture map to be used as the source of reflectance
information, which is then mapped onto the shape based on its UV
parameterization. When no parameters are specified, the model uses the default
of 50% reflectance.

Note that this material is one-sided---that is, observed from the
back side, it will be completely black. If this is undesirable,
consider using the :ref:`twosided <bsdf-twosided>` BRDF adapter plugin.
The following XML snippet describes a diffuse material,
whose reflectance is specified as an sRGB color:

.. code-block:: xml
    :name: diffuse-srgb

    <bsdf type="diffuse">
        <rgb name="reflectance" value="0.2, 0.25, 0.7"/>
    </bsdf>

Alternatively, the reflectance can be textured:

.. code-block:: xml
    :name: diffuse-texture

    <bsdf type="diffuse">
        <texture type="bitmap" name="reflectance">
            <string name="filename" value="wood.jpg"/>
        </texture>
    </bsdf>

*/
template <typename Float, typename Spectrum>
class SmoothDiffuse final : public BSDF<Float, Spectrum> {
public:
    MTS_IMPORT_BASE(BSDF, m_flags, m_components)
    MTS_IMPORT_TYPES(Texture)

    SmoothDiffuse(const Properties &props) : Base(props) {
        m_reflectance = props.texture<Texture>("reflectance", .5f);
        m_flags = BSDFFlags::DiffuseReflection | BSDFFlags::FrontSide;
        m_components.push_back(m_flags);
    }

    MTS_INLINE Spectrum get_spectrum(const SurfaceInteraction3f &si) const {
        Spectrum wavelengths;
        if constexpr (is_spectral_v<Spectrum>) {
            wavelengths[0] = si.wavelengths[0];
            wavelengths[1] = si.wavelengths[1];
            wavelengths[2] = si.wavelengths[2];
            wavelengths[3] = si.wavelengths[3];
        } else
            Throw("Only spectral variants supported");
        return wavelengths;
    }

    size_t foo(const UnpolarizedSpectrum &wavelengths) const {
        auto n_indices = zero<Array<Complex<UnpolarizedSpectrum>, 51>>();
        auto distances = zero<Array<UnpolarizedSpectrum, 51>>();

        UnpolarizedSpectrum test = wavelengths;

        auto k_cortex_indices =
            load<Array<uint32_t, k_cortex_count>>(k_cortex_indices_.data());

        auto k_separation_indices =
            load<Array<uint32_t, k_separation_count>>(k_separation_indices_.data());

        auto m_membrane_indices =
            load<Array<uint32_t, m_membrane_count>>(m_membrane_indices_.data());

        auto m_internal_indices =
            load<Array<uint32_t, m_internal_count>>(m_internal_indices_.data());

        auto m_internal_top_indices =
            load<Array<uint32_t, m_internal_top_count>>(m_internal_top_indices_.data());

        auto air_indices =
            load<Array<uint32_t, air_count>>(air_indices_.data());

        // // Thicknesses in nm
        // Packet<ScalarFloat, k_cortex_count> k_cortex_thickness(5.0f);
        // Packet<ScalarFloat, k_separation_count> k_separation_thickness(50.f);
        // Packet<ScalarFloat, m_membrane_count> m_membrane_thickness(30.f);
        // // Note: this can be changed to adjust periodicity
        // Packet<ScalarFloat, m_internal_count> m_internal_thickness(100.0f);
        // Packet<ScalarFloat, m_internal_top_count> m_internal_top_thickness(50.0f);
        // Packet<ScalarFloat, air_count> air_thickness(std::numeric_limits<ScalarFloat>::max());

        // scatter(&distances, k_cortex_thickness, k_cortex_indices);
        // scatter(&distances, k_separation_thickness, k_separation_indices);
        // scatter(&distances, m_membrane_thickness, m_membrane_indices);
        // scatter(&distances, m_internal_thickness, m_internal_indices);
        // scatter(&distances, m_internal_top_thickness, m_internal_top_indices);
        // scatter(&distances, air_thickness, air_indices);

        // distances[0] = UnpolarizedSpectrum(std::numeric_limits<ScalarFloat>::max());
        // distances[1] = UnpolarizedSpectrum(5.0f);

        // auto x = zero<Array<UnpolarizedSpectrum, 5>>();
        // UInt32 idx = arange<UInt32>() / 2;
        // scatter(&x, full<Array<UnpolarizedSpectrum>>(1.0f), indx);
        // std::cout << idx << '\n';
        // std::cout << x << '\n';

        // using T = Array<UnpolarizedSpectrum>;
        // using Value = value_t<T>;
        // using UInt32P = Packet<uint32_t, T::Size>;

        std::array<uint32_t, Array<UnpolarizedSpectrum>::Size> array = { 9 };
        auto mem = zero<Complex<Array<UnpolarizedSpectrum, 10>>>();
        auto z = sequence<UInt32>(array);
        auto yy = Complex<Array<UnpolarizedSpectrum>>(1.0f, 2.0f);
        scatter(&mem, yy, z);
        std::cout << mem << '\n';
        std::cout << z << '\n';
        std::cout << yy << '\n';

        // auto mem = zero<Array<UnpolarizedSpectrum, 5>>();
        // UInt32 idx(1);
        // scatter(&mem, full<Array<UnpolarizedSpectrum>>(1.0f), idx);
        // std::cout << idx << '\n';
        // std::cout << mem << '\n';
        // auto y = UInt32::Size;

        std::cout << k_cortex_indices << '\n';
        std::cout << test << '\n';
        std::cout << distances << '\n';

        return 5;
    }

    std::pair<BSDFSample3f, Spectrum> sample(const BSDFContext &ctx,
                                             const SurfaceInteraction3f &si,
                                             Float /* sample1 */,
                                             const Point2f &sample2,
                                             Mask active) const override {
        MTS_MASKED_FUNCTION(ProfilerPhase::BSDFSample, active);

        Float cos_theta_i = Frame3f::cos_theta(si.wi);
        BSDFSample3f bs = zero<BSDFSample3f>();

        active &= cos_theta_i > 0.f;
        if (unlikely(none_or<false>(active) ||
                     !ctx.is_enabled(BSDFFlags::DiffuseReflection)))
            return { bs, 0.f };

        bs.wo = warp::square_to_cosine_hemisphere(sample2);
        bs.pdf = warp::square_to_cosine_hemisphere_pdf(bs.wo);
        bs.eta = 1.f;
        bs.sampled_type = +BSDFFlags::DiffuseReflection;
        bs.sampled_component = 0;

        UnpolarizedSpectrum value = m_reflectance->eval(si, active);

        return { bs, select(active && bs.pdf > 0.f, unpolarized<Spectrum>(value), 0.f) };
    }

    Spectrum eval(const BSDFContext &ctx, const SurfaceInteraction3f &si,
                  const Vector3f &wo, Mask active) const override {
        MTS_MASKED_FUNCTION(ProfilerPhase::BSDFEvaluate, active);

        if (!ctx.is_enabled(BSDFFlags::DiffuseReflection))
            return 0.f;

        Float cos_theta_i = Frame3f::cos_theta(si.wi),
              cos_theta_o = Frame3f::cos_theta(wo);

        active &= cos_theta_i > 0.f && cos_theta_o > 0.f;

        UnpolarizedSpectrum value =
            m_reflectance->eval(si, active) * math::InvPi<Float> * cos_theta_o;

        // enum Color { red, green, blue };
        // using MyFloat = mask_t<Color>;
        // using FloatP = Array<float, 3>;
        // FloatP x(1, 2, 3);
        // Vector4f y = arange<Vector4f>() + 3;

        // y[eq(x, Color::red)] += 5.0f;

        // Log(Info, "Value of index 0 : %d", y[0]);

        // auto z = Eigen::Infinity;

        auto x = foo(unpolarized<Spectrum>(get_spectrum(si)));

        return select(active, unpolarized<Spectrum>(value), 0.f);
    }

    Float pdf(const BSDFContext &ctx, const SurfaceInteraction3f &si,
              const Vector3f &wo, Mask active) const override {
        MTS_MASKED_FUNCTION(ProfilerPhase::BSDFEvaluate, active);

        if (!ctx.is_enabled(BSDFFlags::DiffuseReflection))
            return 0.f;

        Float cos_theta_i = Frame3f::cos_theta(si.wi),
              cos_theta_o = Frame3f::cos_theta(wo);

        Float pdf = warp::square_to_cosine_hemisphere_pdf(wo);

        return select(cos_theta_i > 0.f && cos_theta_o > 0.f, pdf, 0.f);
    }

    void traverse(TraversalCallback *callback) override {
        callback->put_object("reflectance", m_reflectance.get());
    }

    std::string to_string() const override {
        std::ostringstream oss;
        oss << "SmoothDiffuse[" << std::endl
            << "  reflectance = " << string::indent(m_reflectance) << std::endl
            << "]";
        return oss.str();
    }

    MTS_DECLARE_CLASS()
private:
    ref<Texture> m_reflectance;
};

MTS_IMPLEMENT_CLASS_VARIANT(SmoothDiffuse, BSDF)
MTS_EXPORT_PLUGIN(SmoothDiffuse, "Smooth diffuse material")
NAMESPACE_END(mitsuba)
