/**
 * @file progress_bar.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implements a progress bar
 * @version 1.2
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

#include "progress_bar.hpp"

namespace CLONES
{

namespace UI
{

ProgressBar::ProgressBar():
    ProgressBar(std::cout, true)
{}

ProgressBar::ProgressBar(std::ostream& out):
    ProgressBar(out, false)
{}

#if WITH_INDICATORS
void ProgressBar::build_new_indicator()
{
    using namespace indicators;

    if (indicator != nullptr) {
        delete indicator;
    }

    indicator = new indicators::ProgressBar{
        option::BarWidth{40},
        option::Start{" ["},
        option::Fill{"█"},
        option::Lead{"█"},
        option::Remainder{"-"},
        option::End{"]"},
        option::ForegroundColor{Color::yellow},
        option::ShowElapsedTime{true},
        option::ShowPercentage{true},
        option::FontStyles{std::vector<FontStyle>{FontStyle::bold}},
        option::Stream{output_stream}
    };

    hide_console_cursor();
}
#endif // WITH_INDICATORS

ProgressBar::ProgressBar(std::ostream& out, const bool hide):
    output_stream(out), last_update(std::chrono::system_clock::now()), percentage(0),
    message(), updated(false)
#if WITH_INDICATORS
    , indicator(nullptr)
#endif // WITH_INDICATORS
{
    using namespace std::chrono_literals;

    update_interval = 1s;

    if (!hide) {
#if WITH_INDICATORS
        build_new_indicator();
#endif  // WITH_INDICATORS
    }
}

ProgressBar& ProgressBar::set_message(const std::string& message)
{
    this->message = message;

#if WITH_INDICATORS
    using namespace indicators;

    if (indicator!=nullptr && (percentage != 100 || !updated)) {
        // the following line is a temporary fix for indicators
        // issue #126 (https://github.com/p-ranav/indicators/issues/126)
        std::string curLocale = setlocale(LC_ALL, nullptr);

        indicator->set_option(option::PostfixText{message});
        indicator->set_progress(percentage);

        // the following line is a temporary fix for indicators
        // issue #126 (https://github.com/p-ranav/indicators/issues/126)
        setlocale(LC_ALL, curLocale.c_str());
    }
#endif  // WITH_INDICATORS

    updated = true;

    return *this;
}

ProgressBar& ProgressBar::set_progress(const unsigned int percentage, const std::string& message)
{
    using namespace std::chrono;

    this->percentage = std::min(percentage, static_cast<unsigned int>(100));
    this->message = message;

    updated = false;

    const auto curr_time = system_clock::now();
    if (duration_cast<seconds>(curr_time-last_update).count()>=1) {
        last_update = curr_time;

        set_message(message);
    }

    return *this;
}


ProgressBar& ProgressBar::set_progress(const unsigned int percentage)
{
    using namespace std::chrono;

    this->percentage = percentage;

    updated = false;

    const auto curr_time = system_clock::now();
    if (duration_cast<seconds>(curr_time-last_update).count()>=1) {
        last_update = curr_time;

        updated = true;

#if WITH_INDICATORS
        if (indicator!=nullptr) {
            indicator->set_progress(percentage);
        }
#endif  // WITH_INDICATORS
    }

    return *this;
}

const unsigned int& ProgressBar::get_progress() const
{
    return percentage;
}

void ProgressBar::show_console_cursor()
{
#if WITH_INDICATORS
    indicators::show_console_cursor(true);
#endif  // WITH_INDICATORS
}

void ProgressBar::hide_console_cursor()
{
#if WITH_INDICATORS
    indicators::show_console_cursor(false);
#endif  // WITH_INDICATORS
}

void ProgressBar::close()
{
    set_message(this->message);

#if WITH_INDICATORS
    if (indicator != nullptr) {
        if (percentage != 100) {
            indicator->mark_as_completed();
        }

        delete indicator;

        indicator = nullptr;
    }
#endif  // WITH_INDICATORS

    show_console_cursor();
}

ProgressBar& ProgressBar::init_new()
{

#if WITH_INDICATORS
    if (indicator != nullptr) {
        close();

        build_new_indicator();
    }
#endif  // WITH_INDICATORS

    percentage = 0;
    last_update = std::chrono::system_clock::now();
    message = "";
    updated = false;

    return *this;
}

ProgressBar::~ProgressBar()
{
    close();
}

}  // UI

} // CLONES