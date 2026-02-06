# -*- coding: utf-8 -*-
#
# Copyright \(c\) 2023-2026 Alberto Casagrande <alberto.casagrande@uniud.it>
#
# MIT License
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import unittest

import CLONES


class TestCellEventType(unittest.TestCase):

    def test_values(self):
        self.assertEqual(CLONES.CellEventType.values, {
            0: CLONES.CellEventType.DEATH,
            1: CLONES.CellEventType.DUPLICATION,
            2: CLONES.CellEventType.EPIGENETIC_SWITCH,
            3: CLONES.CellEventType.MUTATION
        })


class TestEpigeneticRates(unittest.TestCase):
    def test_init(self):
        try:
            CLONES.EpigeneticRates(0.2, 0.1)
        except BaseException:
            self.fail("CLONES.EpigeneticRates(0.2, 0.1) unexpectedly"
                      + " raised an exeception!")

        with self.assertRaises(Exception):
            CLONES.EpigeneticRates()

        with self.assertRaises(Exception):
            CLONES.EpigeneticRates(-0.2, 0.1)

        with self.assertRaises(Exception):
            CLONES.EpigeneticRates(0.2, -0.1)

        with self.assertRaises(Exception):
            CLONES.EpigeneticRates(1.2, 0.1)

        with self.assertRaises(Exception):
            CLONES.EpigeneticRates(0.2, 1.1)

    def test_get(self):
        e_rates = CLONES.EpigeneticRates(0.2, 0.1)

        self.assertEqual(e_rates.get_methylation_rate(), 0.2)
        self.assertEqual(e_rates.get_demethylation_rate(), 0.1)

    def test_set(self):
        e_rates = CLONES.EpigeneticRates(0.2, 0.1)

        e_rates.set_methylation_rate(0.3)
        e_rates.set_demethylation_rate(0.4)

        self.assertEqual(e_rates.get_methylation_rate(), 0.3)
        self.assertEqual(e_rates.get_demethylation_rate(), 0.4)

        with self.assertRaises(Exception):
            e_rates.set_methylation_rate(1.3)

        with self.assertRaises(Exception):
            e_rates.set_demethylation_rate(1.3)

        with self.assertRaises(Exception):
            e_rates.set_methylation_rate(-0.3)

        with self.assertRaises(Exception):
            e_rates.set_demethylation_rate(-0.3)


class TestClone(unittest.TestCase):
    def test_init(self):
        try:
            CLONES.Clone("A", [CLONES.EpigeneticRates(0.01, 0.01)])
        except BaseException:
            self.fail('CLONES.Clone("A", '
                      + '[CLONES.EpigeneticRates(0.01, 0.01)]) '
                      + ' raised an unexpected exeception!')

        try:
            CLONES.Clone("A", [[0.01, 0.01]])
        except BaseException:
            self.fail('CLONES.Clone("A", [[0.01, 0.01]]) '
                      + ' raised an unexpected exeception!')

        try:
            CLONES.Clone("A", [])
        except BaseException:
            self.fail('CLONES.Clone("A", []) '
                      + ' raised an unexpected exeception!')

        try:
            CLONES.Clone("A", [[0.01, 0.01], [0.01, 0.01]])
        except BaseException:
            self.fail('CLONES.Clone("A", [[0.01, 0.01], '
                      + '[0.01, 0.01]]) raised an unexpected exeception!')

        with self.assertRaises(Exception):
            CLONES.Clone("A", 2)

        with self.assertRaises(Exception):
            CLONES.Clone("A", ['a'])

        with self.assertRaises(Exception):
            CLONES.Clone("A", [[0.01, 0.01, 2]])

    def test_properties(self):
        A = CLONES.Clone("A", [[0.01, 0.01]])

        self.assertEqual(A.name, "A")
        self.assertEqual(A.num_of_promoters, 1)

        B = CLONES.Clone("B", [])

        self.assertEqual(B.name, "B")
        self.assertEqual(B.id, A.id+1)
        self.assertEqual(B.num_of_promoters, 0)

        with self.assertRaises(Exception):
            A.name = "test"

        with self.assertRaises(Exception):
            A.id = 3

        with self.assertRaises(Exception):
            A.num_of_promoters = 7

    def test_set_rates(self):
        A = CLONES.Clone("A", [[0.01, 0.01]])
        try:
            A.set_rates("-", {CLONES.CellEventType.DEATH: 0.1,
                              CLONES.CellEventType.DUPLICATION: 0.2})
        except BaseException:
            self.fail('A.set_rates("-", {CLONES.CellEventType.DEATH: 0.1, '
                      + 'CLONES.CellEventType.DUPLICATION: 0.2})'
                      + ' raised an unexpected exeception!')

        with self.assertRaises(Exception):
            A.set_rates("", {CLONES.CellEventType.DEATH: 0.1,
                             CLONES.CellEventType.DUPLICATION: 0.2})

        B = CLONES.Clone("B", [])

        try:
            B.set_rates("", {CLONES.CellEventType.DEATH: 0.1,
                             CLONES.CellEventType.DUPLICATION: 0.2})
        except BaseException:
            self.fail('B.set_rates("", {CLONES.CellEventType.DEATH: 0.1, '
                      + 'CLONES.CellEventType.DUPLICATION: 0.2})'
                      + ' raised an unexpected exeception!')

        with self.assertRaises(Exception):
            B.set_rates("-", {CLONES.CellEventType.DEATH: 0.1,
                              CLONES.CellEventType.DUPLICATION: 0.2})

    def test_get_rate(self):
        A = CLONES.Clone("A", [[0.01, 0.01]])
        A.set_rates("-", {CLONES.CellEventType.DEATH: 0.1,
                          CLONES.CellEventType.DUPLICATION: 0.2})
        try:
            self.assertEqual(A.get_rate("-", CLONES.CellEventType.DEATH), 0.1)
            self.assertEqual(A.get_rate("-", CLONES.CellEventType.DUPLICATION),
                             0.2)
            self.assertEqual(A.get_rate("-",
                                        CLONES.CellEventType.EPIGENETIC_SWITCH),
                             0.0)
        except BaseException:
            self.fail('CLONES.Clone.set_rate() '
                      + ' raised an unexpected exeception!')

        with self.assertRaises(Exception):
            self.assertEqual(A.get_rate("", CLONES.CellEventType.DEATH), 0.01)

        B = CLONES.Clone("B", [])
        B.set_rates("", {CLONES.CellEventType.DEATH: 0.1,
                         CLONES.CellEventType.DUPLICATION: 0.2})
        try:
            self.assertEqual(B.get_rate("", CLONES.CellEventType.DEATH), 0.1)
            self.assertEqual(B.get_rate("", CLONES.CellEventType.DUPLICATION),
                             0.2)
            self.assertEqual(B.get_rate("",
                                        CLONES.CellEventType.EPIGENETIC_SWITCH),
                             0.0)
        except BaseException:
            self.fail('CLONES.Clone.set_rate() '
                      + ' raised an unexpected exeception!')

        with self.assertRaises(Exception):
            self.assertEqual(B.get_rate("-", CLONES.CellEventType.DEATH), 0.01)


class TestSimulation(unittest.TestCase):
    def test_init(self):
        try:
            CLONES.Simulation()
        except BaseException:
            self.fail('CLONES.Simulation() raised an unexpected exeception!')

        try:
            CLONES.Simulation(5)
        except BaseException:
            self.fail('CLONES.Simulation(5) raised an unexpected exeception!')

        try:
            CLONES.Simulation(5, 1)
        except BaseException:
            self.fail('CLONES.Simulation(5, 1) raised an unexpected'
                      + ' exeception!')

    def test_set_tissue(self):
        sim = CLONES.Simulation()

        try:
            sim.set_tissue("Liver", [100, 100])
        except BaseException:
            self.fail('sim.set_tissue("Liver", [100, 100]) raised'
                      + ' an unexpected exeception!')

        try:
            sim.set_tissue("Liver", [100, 100, 100])
        except BaseException:
            self.fail('sim.set_tissue("Liver", [100, 100, 100]) raised'
                      + ' an unexpected exeception!')

        with self.assertRaises(Exception):
            sim.set_tissue("Liver", [100])

        with self.assertRaises(Exception):
            sim.set_tissue("Liver", [])

        with self.assertRaises(Exception):
            sim.set_tissue("Liver", [100, 1, 1, 1])

        with self.assertRaises(Exception):
            sim.set_tissue("Liver", [100, 100, -100])

        with self.assertRaises(Exception):
            sim.set_tissue("Liver", [100, -100])

    def test_get_time(self):
        sim = CLONES.Simulation()

        self.assertEqual(sim.get_time(), 0)

    def test_add_mutant(self):
        sim = CLONES.Simulation()

        A = CLONES.Clone("A", [[0.01, 0.01]])
        A.set_rates("-", {CLONES.CellEventType.DEATH: 0.1,
                          CLONES.CellEventType.DUPLICATION: 0.2})
        A.set_rates("+", {CLONES.CellEventType.DEATH: 0.01,
                          CLONES.CellEventType.DUPLICATION: 0.02})

        sim.set_tissue("Liver", [100, 100])

        try:
            # now the simulation has a tissue
            sim.add_mutant(A)

        except BaseException:
            self.fail('sim.add_mutant(A) raised an unexpected exeception!')

        with self.assertRaises(Exception):
            sim.add_mutant("A")

    def test_place_cell(self):
        sim = CLONES.Simulation()

        A = CLONES.Clone("A", [[0.01, 0.01]])
        A.set_rates("-", {CLONES.CellEventType.DEATH: 0.1,
                          CLONES.CellEventType.DUPLICATION: 0.2})
        A.set_rates("+", {CLONES.CellEventType.DEATH: 0.01,
                          CLONES.CellEventType.DUPLICATION: 0.02})

        sim.set_tissue("Liver", [100, 100])

        sim.add_mutant(A)

        try:
            sim.place_cell(A, "-", [50, 50])

        except BaseException:
            self.fail('sim.place_cell(A, "-", [50, 50]) raised ' +
                      'an unexpected exeception!')

        with self.assertRaises(Exception):
            sim.place_cell(A, "-", [50, 150])

    def test_schedule_mutation(self):
        sim = CLONES.Simulation()

        sim.set_tissue("Liver", [100, 100])

        A = CLONES.Clone("A", [[0.01, 0.01]])
        A.set_rates("-", {CLONES.CellEventType.DEATH: 0.1,
                          CLONES.CellEventType.DUPLICATION: 0.2})
        A.set_rates("+", {CLONES.CellEventType.DEATH: 0.01,
                          CLONES.CellEventType.DUPLICATION: 0.02})
        sim.add_mutant(A)

        B = CLONES.Clone("B", [[0.01, 0.01]])
        B.set_rates("-", {CLONES.CellEventType.DEATH: 0.1,
                          CLONES.CellEventType.DUPLICATION: 0.3})
        B.set_rates("+", {CLONES.CellEventType.DEATH: 0.1,
                          CLONES.CellEventType.DUPLICATION: 0.45})
        sim.add_mutant(B)

        C = CLONES.Clone("C", [])
        C.set_rates("", {CLONES.CellEventType.DEATH: 0.1,
                         CLONES.CellEventType.DUPLICATION: 0.3})
        sim.add_mutant(C)

        try:
            sim.schedule_mutation(A, B, 70)
        except BaseException:
            self.fail('sim.schedule_mutation(A, B, 70) raised'
                      + ' an unexpected exeception!')

        with self.assertRaises(Exception):
            # methylation signature incompatible
            sim.schedule_mutation(A, C, 70)

        try:
            # methylation signature compatible
            sim.schedule_mutation(C, B, 70)
        except BaseException:
            self.fail('sim.schedule_mutation(C, B, 70) raised'
                      + ' an unexpected exeception!')

    def test_run_up_to(self):
        sim = CLONES.Simulation()

        sim.set_tissue("Liver", [100, 100])

        A = CLONES.Clone("A", [[0.01, 0.01]])
        A.set_rates("-", {CLONES.CellEventType.DEATH: 0.1,
                          CLONES.CellEventType.DUPLICATION: 0.2})
        A.set_rates("+", {CLONES.CellEventType.DEATH: 0.01,
                          CLONES.CellEventType.DUPLICATION: 0.02})
        sim.add_mutant(A)

        B = CLONES.Clone("B", [[0.01, 0.01]])
        B.set_rates("-", {CLONES.CellEventType.DEATH: 0.1,
                          CLONES.CellEventType.DUPLICATION: 0.3})
        B.set_rates("+", {CLONES.CellEventType.DEATH: 0.1,
                          CLONES.CellEventType.DUPLICATION: 0.45})
        sim.add_mutant(B)
        sim.schedule_mutation(A, B, 70)

        sim.storage_enabled = False

        sim.run_up_to(0, quiet=True)

    def test_death_activation_level(self):
        sim = CLONES.Simulation()

        try:
            value = sim.death_activation_level
        except BaseException:
            self.fail('reading sim.death_activation_level raised'
                      + ' an unexpected exeception!')

        self.assertTrue(isinstance(value, int))

        try:
            sim.death_activation_level = 100
        except BaseException:
            self.fail('setting sim.death_activation_level raised'
                      + ' an unexpected exeception!')

        self.assertEqual(sim.death_activation_level, 100)


if __name__ == '__main__':
    unittest.main()
