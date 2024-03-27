
        // for (size_t i = 0; i < 2; ++i)
        // {
        //     std::cout << airIndices[i] << '\n';
        // }

        std::array<uint32_t, 2> test_indices = {1, 2};
        static constexpr size_t test_num_layers = std::size(test_indices);

        Array<uint32_t, test_num_layers> testIndices =
            load<Array<uint32_t, test_num_layers>>(test_indices.data());

        Array<ScalarFloat, test_num_layers> test_thickness(
            std::numeric_limits<float>::max());

        scatter(&distances, test_thickness, testIndices);

        // Array<float, 3> arr(1.f);
        // Array<uint32_t, 3> idx(1, 2, 3);
        // Array<UnpolarizedSpectrum, 4> dest = zero<Array<UnpolarizedSpectrum, 4>>();
        // scatter(&dest, arr, idx);