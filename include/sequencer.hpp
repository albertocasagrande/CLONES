/**
 * @file sequencer.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines sequencer models
 * @version 1.8
 * @date 2026-02-06
 *
 * @copyright Copyright (c) 2023-2026
 *
 * MIT License
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#ifndef __CLONES_SEQUENCER__
#define __CLONES_SEQUENCER__

#include <string>
#include <random>
#include <vector>
#include <cmath>

#include "genomic_position.hpp"
#include "genomic_sequence.hpp"
#include "read.hpp"

namespace CLONES
{

/**
 * @brief The namespace of sequencer models
 */
namespace Sequencers
{

/**
 * @brief A basic sequencer model
 */
class BasicSequencer
{
public:

    /**
     * @brief Get the model name
     *
     * @return the model name
     */
    virtual std::string get_model_name() const = 0;

    /**
     * @brief Get the platform name
     *
     * This method returns the sequencer platform as specified by [1].
     *
     * [1] "Sequence Alignment/Map Format Specification", The SAM/BAM Format Specification
     *     Working Group, 22 August 2022, http://samtools.github.io/hts-specs/SAMv1.pdf
     *
     * @return the platform name
     */
    virtual std::string get_platform_name() const = 0;

    /**
     * @brief Simulate the sequencing of a base
     *
     *
     * This method simulates the sequencing of a base. It applies an error to the
     * base according to the sequencer model and it also produces the corresponding
     * quality scores.
     *
     * @param[in,out] base is the base whose sequencing must be simulated
     * @param[in] read_size is the simulated read size
     * @param[in] base_position is the position of the base in the read
     * @param[in] read_position is the position of the first read base in the reference
     *      genome
     * @param[in] reverse is a Boolean flag for simulating reverse read
     * @return the quality score of the base in Sanger FASTQ format
     */
    virtual char simulate_seq(char& base, const size_t& read_size,
                              const size_t base_position,
                              const Mutations::GenomicPosition& read_position,
                              const bool& reverse=false) = 0;

    /**
     * @brief Simulate the sequencing of a read
     *
     * This method simulates the sequencing of a read. It applies some errors to the
     * read according to the sequencer model and it also produces the corresponding
     * quality scores.
     *
     * @param[in,out] read is the read whose sequencing must be simulated
     * @param[in] position is the position of the first read base in the reference
     *      genome
     * @param[in] reverse is a Boolean flag for simulating reverse read
     * @return the quality scores of the read in Sanger FASTQ format
     */
    virtual std::string simulate_seq(Mutations::SequencingSimulations::Read& read,
                                     const Mutations::GenomicPosition& position,
                                     const bool& reverse=false) = 0;

    /**
     * @brief Test whether paired reads are supported by the sequencer
     *
     * @return `true` if and only if the sequencer supports paired reads
     */
    virtual inline bool supports_paired_reads() const
    {
        return false;
    }

    /**
     * @brief The destroyer
     */
    virtual inline ~BasicSequencer()
    {}
};

/**
 * @brief The Sanger quality codec
 *
 * This class implements the methods to encode and decode a base
 * error probability in a from a quality code as specified by
 * Sanger standard.
 */
class SangerQualityCodec
{
public:
    /**
     * @brief Get the base quality code
     *
     * @return The base quality code
     */
    inline constexpr static uint8_t base_quality_code()
    {
        return 33;
    }

    /**
     * @brief Get the minimum quality value
     *
     * @return The minimum quality value
     */
    inline constexpr static uint8_t min_quality_value()
    {
        return 0;
    }

    /**
     * @brief Get the maximum quality value
     *
     * @return The maximum quality value
     */
    inline constexpr static uint8_t max_quality_value()
    {
        // the maximum quality value in phred is 40, however,
        // the quality value never tops 37
        return 37;
    }

    /**
     * @brief Get the minimum quality code
     *
     * @return The minimum quality code
     */
    inline constexpr static uint8_t min_quality_code()
    {
        return min_quality_value() + base_quality_code();
    }

    /**
     * @brief Get the maximum quality code
     *
     * @return The maximum quality code
     */
    inline constexpr static uint8_t max_quality_code()
    {
        return max_quality_value() + base_quality_code();
    }

    /**
     * @brief Encode an error probability into a quality value
     *
     * @param error_prob is the error probability to be encoded
     * @return the quality code corresponding to `error_prob`
     */
    static char encode(const double error_prob);

    /**
     * @brief Decode a quality code into an error probability
     *
     * @param quality_code is the quality code to be decoded
     * @return the error probability corresponding to `quality_code`
     */
    static double decode(const char quality_code);
};

/**
 * @brief A statistical quality score model
 *
 * This statistical quality score model simulates read quality scores,
 * exclusively taking into account the base positions in the simulated read.
 * The scores are sampled on normal distributions whose means and standard
 * deviations are estimated by using two functions whose domains are the
 * read positions.
 */

template<typename QUALITY_CODEC=SangerQualityCodec>
class QualityScoreModel
{
public:
    using EstimatingFunctionType = std::function<double(const size_t&)>;

private:
    EstimatingFunctionType mean;   //!< The estimator of the mean quality score
    EstimatingFunctionType stddev; //!< The estimator of the standard deviation quality score

    std::vector<std::normal_distribution<double>> quality_dists; //!< The quality distributions

    /**
     * @brief A function estimating the mean quality score on a read position
     *
     * The mean quality values of Illumina sequencers per position tend to
     * linearly increase from the first base reaching a peak around the 7th
     * base and, then, linearly sink. Because of this, this function models
     * the mean by the linear switching function
     * $$
     * f_m(x) = \begin{cases}
     * a_1 x + b & \textrm{if $x<p$}\\
     * a_2 x + b + p (a_1-a_2) & \textrm{if $x\geq p$}
     * \end{cases}
     * $$
     * where \$p\$ is the mean peak position, \$b\$ is the mean quality
     * score for the first base, and \$a_1\$ and \$a_2\$ are the linear
     * coefficients before and after the peak, respectively.
     * The parameters were fitted to real sequencing data.
     *
     * @param position is the position for which the mean quality score
     *  is requested
     * @return The mean quality score on `position`
     */
    static double default_mean(const size_t& position)
    {
        const double library_const = static_cast<double>(100)/150;
        constexpr double mean_a1 = 0.0845729; //!< The mean linear coefficient before the peak
        const double mean_a2 = -0.0146604*library_const; //!< The mean linear coefficient after the peak
        constexpr double mean_b = 35.1498051;  //!< The mean quality score in first position

        constexpr uint16_t mean_peak_pos = 7; //!< The position of the mean peak

        if (position<mean_peak_pos) {
            return mean_a1 * position + mean_b;
        }

        return mean_a2*position + (mean_peak_pos*(mean_a1-mean_a2)+mean_b);
    }

    /**
     * @brief A function estimating the standard deviation quality score on a read position
     *
     * This function models the standard deviations of the quality scores by the quadratic
     * function
     * $$
     * f_s(x) = c x^2 + d.
     * $$
     * The parameters were fitted to real sequencing data.
     *
     * @param position is the position for which the standard deviation
     *   quality score is requested
     * @return The standard deviation quality score on `position`
     */
    static double default_stddev(const size_t& position)
    {
        const double library_const = static_cast<double>(100)/150;

        const double stddev_quad_coeff = 9.661484e-06 * library_const * library_const;   //!< The standard deviation quadratic coefficient
        constexpr double stddev_base_coeff = 1.731231e-01;   //!< The standard deviation zero-degree coefficient

        return stddev_quad_coeff*position*position + stddev_base_coeff;
    }

public:
    using QualityScoreCodec = QUALITY_CODEC;

    /**
     * @brief Construct a new quality score model
     *
     * @param mean_qual_function is the estimator of the mean quality score in a position
     * @param stddev_qual_function
     */
    explicit QualityScoreModel(const EstimatingFunctionType& mean_qual_function=default_mean,
                               const EstimatingFunctionType& stddev_qual_function=default_stddev):
         mean(mean_qual_function), stddev(stddev_qual_function)
    {
        if constexpr(!std::is_base_of_v<QUALITY_CODEC, SangerQualityCodec>) {
            throw std::logic_error("This class only supports `SangerQualityCodec` at the moment");
        }
    }

    /**
     * @brief Simulate the quality score of a base
     *
     * @tparam RANDOM_GENERATOR is the type of the pseudo-random number generator
     * @param generator is a pseudo-random number generator
     * @param position is the read position of the base for which the quality
     *   score is simulated
     * @return The simulated quality score
     */
    template<typename RANDOM_GENERATOR>
    uint8_t operator()(RANDOM_GENERATOR& generator, const size_t& position)
    {
        while (position >= quality_dists.size()) {
            quality_dists.emplace_back(mean(quality_dists.size()),
                                       stddev(quality_dists.size()));
        }

        const auto quality = quality_dists[position](generator);

        if (quality > QUALITY_CODEC::max_quality_value()) {
            return QUALITY_CODEC::max_quality_code();
        }

        if (quality < QUALITY_CODEC::min_quality_value()) {
            return QUALITY_CODEC::min_quality_code();
        }

        return static_cast<uint8_t>(quality) + QUALITY_CODEC::base_quality_code();
    }
};

/**
 * @brief A constant quality score model
 *
 * This model assigns the same score to any base in any
 * position of a simulated read.
 *
 * @tparam QUALITY_CODEC is the quality score codec
 */

template<typename QUALITY_CODEC=SangerQualityCodec>
class ConstantQualityScoreModel
{
    uint8_t quality_score;  //!< The constant quality score
public:
    using QualityScoreCodec = QUALITY_CODEC;

    /**
     * @brief The empty constructor
     */
    ConstantQualityScoreModel():
        quality_score(QUALITY_CODEC::max_quality_code())
    {}

    /**
     * @brief A constructor
     *
     * @param quality_score is the constant quality score
     */
    explicit ConstantQualityScoreModel(const uint8_t quality_score):
            quality_score(std::max(std::min(quality_score,
                                   QUALITY_CODEC::max_quality_code()),
                          QUALITY_CODEC::min_quality_code()))
    {}

    /**
     * @brief Simulate the quality score of a base
     *
     * @tparam RANDOM_GENERATOR is the type of the pseudo-random number generator
     * @param generator is a pseudo-random number generator
     * @param position is the read position of the base for which the quality
     *   score is simulated
     * @return The constant quality score `ConstantQualityScoreModel::quality_score`
     */
    template<typename RANDOM_GENERATOR>
    inline uint8_t operator()(RANDOM_GENERATOR& generator, const size_t& position)
    {
        (void)generator;
        (void)position;

        return quality_score;
    }
};

/**
 * @brief A namespace for Illumina sequencer
 */
namespace Illumina
{

/**
 * @brief A basic error-less Illumina sequencer
 *
 * @tparam QUALITY_CODEC is the type of the quality score codec
 */
template<typename QUALITY_CODEC = SangerQualityCodec>
class ErrorLessSequencer : public CLONES::Sequencers::BasicSequencer
{
public:
    using QualityScoreCodec = QUALITY_CODEC;

    /**
     * @brief Construct a new error-less Illumina sequencer object
     */
    ErrorLessSequencer()
    {}

    /**
     * @brief Get the model name
     *
     * @return the model name
     */
    inline std::string get_model_name() const override
    {
        return "Error-less Illumina Sequencer";
    }

    /**
     * @brief Get the platform name
     *
     * This method returns the sequencer platform as specified by [1].
     *
     * [1] "Sequence Alignment/Map Format Specification", The SAM/BAM Format Specification
     *     Working Group, 22 August 2022, http://samtools.github.io/hts-specs/SAMv1.pdf
     *
     * @return the platform name
     */
    inline std::string get_platform_name() const override
    {
        return "ILLUMINA";
    }

    /**
     * @brief Simulate the sequencing of a base
     *
     *
     * This method simulates the sequencing of a base. It applies an error to the
     * base according to the sequencer model and it also produces the corresponding
     * quality scores.
     *
     * @param[in,out] base is the base whose sequencing must be simulated
     * @param[in] read_size is the simulated read size
     * @param[in] base_position is the position of the base in the read
     * @param[in] read_position is the position of the first read base in the reference
     *      genome
     * @param[in] reverse is a Boolean flag for simulating reverse read
     * @return the quality score of the base in Sanger FASTQ format
     */
    char simulate_seq(char& base, const size_t& read_size,
                      const size_t base_position,
                      const Mutations::GenomicPosition& read_position,
                      const bool& reverse=false) override
    {
        (void)read_size;
        (void)base_position;
        (void)read_position;
        (void)reverse;

        if (base != 'N') {
            return QUALITY_CODEC::max_quality_code();
        }

        return QUALITY_CODEC::min_quality_code();
    }

    /**
     * @brief Simulate the sequencing of a read
     *
     * This method simulates the sequencing of a read. It applies some errors to the
     * read according to the sequencer model and it also produces the corresponding
     * quality scores.
     *
     * @param[in,out] read is the read whose sequencing must be simulated
     * @param[in] position is the position of the first read base in the reference
     *      genome
     * @param[in] reverse is a Boolean flag for simulating reverse read
     * @return the quality scores of the read in Sanger FASTQ format
     */
    std::string simulate_seq(Mutations::SequencingSimulations::Read& read,
                             const Mutations::GenomicPosition& position,
                             const bool& reverse=false) override
    {
        (void)position;
        (void)reverse;

        std::string qual(read.size(), QUALITY_CODEC::max_quality_code());
        for (size_t i=0; i<read.size(); ++i) {
            if (read[i] == 'N') {
                qual[i] = QUALITY_CODEC::min_quality_code();
            }
        }

        return qual;
    }

    /**
     * @brief Test whether paired reads are supported by the sequencer
     *
     * @return `true` if and only if the sequencer supports paired reads
     */
    inline bool supports_paired_reads() const override
    {
        return true;
    }

    /**
     * @brief The destroyer
     */
    ~ErrorLessSequencer()
    {}
};

/**
 * @brief A basic implementation of Illumina sequencer model
 *
 * @tparam QUALITY_SCORE_MODEL is the type of the quality score model
 * @tparam QUALITY_CODEC is the type of the quality value codec
 * @tparam RANDOM_GENERATOR is the type of the random generator
 * This class simulates the sequencing errors by using an error probability
 * that is independent from both the base position in the read and the
 * position in the reference genome. The method uniformly samples a value
 * in the interval [0,1]; whenever the sampled value is smaller than
 * the specified probability threshold a sequencing error is produced.
 * The quality score of a base is determined by using the sample value.
 */
template<template<class> typename QUALITY_SCORE_MODEL = QualityScoreModel,
         typename QUALITY_CODEC = SangerQualityCodec,
         typename RANDOM_GENERATOR = std::mt19937_64>
class BasicSequencer : public ErrorLessSequencer<QUALITY_CODEC>
{
    double error_rate;  //!< The error rate

    QUALITY_SCORE_MODEL<QUALITY_CODEC> quality_score_model;    //!< The quality score model

    RANDOM_GENERATOR number_generator;          //!< The pseudo-random number generator

    std::uniform_real_distribution<double> error_dist;  //!< The uniform error distribution

    /**
     * @brief Get a wrong base
     *
     * @param orig is the original base
     * @return an uniformly selected base different
     *      from `orig`
     */
    char get_wrong_base(const char& orig)
    {
        std::uniform_int_distribution<size_t> offset_dist(1,3);

        auto pos = (GenomicSequence::get_base_index(orig)
                    + offset_dist(number_generator))%4;

        return GenomicSequence::DNA_bases[pos];
    }
public:
    using QualityScoreCodec = QUALITY_CODEC;

    /**
     * @brief Construct a new basic Illumina sequencer object
     *
     * @param error_rate is the error probability in each base
     * @param quality_score_model is the quality score model
     * @param seed is the random number generator seed
     */
    explicit BasicSequencer(const double& error_rate,
                            const QUALITY_SCORE_MODEL<QUALITY_CODEC>& quality_score_model,
                            const int& seed=0):
        error_rate(error_rate), quality_score_model(quality_score_model),
        number_generator(seed), error_dist(0,1)
    {}

    /**
     * @brief Construct a new basic Illumina sequencer object
     *
     * @param error_rate is the error probability in each base
     * @param seed is the random number generator seed
     */
    explicit BasicSequencer(const double& error_rate,
                            const int& seed=0):
        error_rate(error_rate), number_generator(seed), error_dist(0,1)
    {}

    /**
     * @brief Get the model name
     *
     * @return the model name
     */
    inline std::string get_model_name() const override
    {
        return "Basic Illumina Sequencer";
    }

    /**
     * @brief Simulate the sequencing of a base
     *
     *
     * This method simulates the sequencing of a base. It applies an error to the
     * base according to the sequencer model and it also produces the corresponding
     * quality scores.
     *
     * @param[in,out] base is the base whose sequencing must be simulated
     * @param[in] read_size is the simulated read size
     * @param[in] base_position is the position of the base in the read
     * @param[in] read_position is the position of the first read base in the reference
     *      genome
     * @param[in] reverse is a Boolean flag for simulating reverse read
     * @return the quality score of the base in Sanger FASTQ format
     */
    char simulate_seq(char& base, const size_t& read_size,
                      const size_t base_position,
                      const Mutations::GenomicPosition& read_position,
                      const bool& reverse=false) override
    {
        (void)base_position;
        (void)read_position;

        if (base == 'N') {
            return QUALITY_CODEC::min_quality_code();
        }

        auto error_sample = error_dist(number_generator);

        if (error_sample > error_rate) {
            const size_t seq_pos = (reverse?read_size-base_position-1:
                                            base_position);
            return quality_score_model(number_generator, seq_pos);
        }

        base = get_wrong_base(base);

        return QUALITY_CODEC::encode(error_sample/error_rate);
    }

    /**
     * @brief Simulate the sequencing of a read
     *
     * This method simulates the sequencing of a read. It applies some errors to the
     * read according to the sequencer model and it also produces the corresponding
     * quality scores.
     *
     * @param[in,out] read is the read whose sequencing must be simulated
     * @param[in] position is the position of the first read base in the reference
     *      genome
     * @return the quality scores of the read in Sanger FASTQ format
     */
    std::string simulate_seq(Mutations::SequencingSimulations::Read& read,
                             const Mutations::GenomicPosition& position,
                             const bool& reverse=false) override
    {
        std::string qual(read.size(), QUALITY_CODEC::max_quality_code());

        const std::string& read_seq = read.get_sequence();
        auto qual_it = qual.begin();
        auto read_it = read_seq.begin();
        for (size_t base_pos = 0; base_pos < read_seq.size();
                ++base_pos,++qual_it,++read_it) {
            char base = *read_it;
            *qual_it = simulate_seq(base, read_seq.size(), base_pos,
                                    position, reverse);

            read.alter_base(base_pos, base);
        }

        return qual;
    }

    /**
     * @brief Get the error probability per base
     *
     * @return a constant reference to the error probability per base
     */
    inline const double& get_error_rate() const
    {
        return error_rate;
    }

    /**
     * @brief The destroyer
     */
    ~BasicSequencer()
    {}
};

}   // Illumina

}   // Sequencers

}   // CLONES

#endif // __CLONES_SEQUENCER__