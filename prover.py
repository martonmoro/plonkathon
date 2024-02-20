from compiler.program import Program, CommonPreprocessedInput
from utils import *
from setup import *
from typing import Optional
from dataclasses import dataclass
from transcript import Transcript, Message1, Message2, Message3, Message4, Message5
from poly import Polynomial, Basis


@dataclass
class Proof:
    msg_1: Message1
    msg_2: Message2
    msg_3: Message3
    msg_4: Message4
    msg_5: Message5

    def flatten(self):
        proof = {}
        proof["a_1"] = self.msg_1.a_1
        proof["b_1"] = self.msg_1.b_1
        proof["c_1"] = self.msg_1.c_1
        proof["z_1"] = self.msg_2.z_1
        proof["t_lo_1"] = self.msg_3.t_lo_1
        proof["t_mid_1"] = self.msg_3.t_mid_1
        proof["t_hi_1"] = self.msg_3.t_hi_1
        proof["a_eval"] = self.msg_4.a_eval
        proof["b_eval"] = self.msg_4.b_eval
        proof["c_eval"] = self.msg_4.c_eval
        proof["s1_eval"] = self.msg_4.s1_eval
        proof["s2_eval"] = self.msg_4.s2_eval
        proof["z_shifted_eval"] = self.msg_4.z_shifted_eval
        proof["W_z_1"] = self.msg_5.W_z_1
        proof["W_zw_1"] = self.msg_5.W_zw_1
        return proof


@dataclass
class Prover:
    group_order: int
    setup: Setup
    program: Program
    pk: CommonPreprocessedInput

    def __init__(self, setup: Setup, program: Program):
        self.group_order = program.group_order
        self.setup = setup
        self.program = program
        self.pk = program.common_preprocessed_input()

    def prove(self, witness: dict[Optional[str], int]) -> Proof:
        # Initialise Fiat-Shamir transcript
        transcript = Transcript(b"plonk")

        # Collect fixed and public information
        # FIXME: Hash pk and PI into transcript
        public_vars = self.program.get_public_assignments()
        PI = Polynomial(
            [Scalar(-witness[v]) for v in public_vars]
            + [Scalar(0) for _ in range(self.group_order - len(public_vars))],
            Basis.LAGRANGE,
        )
        self.PI = PI

        # Round 1
        msg_1 = self.round_1(witness)
        self.beta, self.gamma = transcript.round_1(msg_1)

        # Round 2
        msg_2 = self.round_2()
        self.alpha, self.fft_cofactor = transcript.round_2(msg_2)

        # Round 3
        msg_3 = self.round_3()
        self.zeta = transcript.round_3(msg_3)

        # Round 4
        msg_4 = self.round_4()
        self.v = transcript.round_4(msg_4)

        # Round 5
        msg_5 = self.round_5()

        return Proof(msg_1, msg_2, msg_3, msg_4, msg_5)

    def round_1(
        self,
        witness: dict[Optional[str], int],
    ) -> Message1:
        program = self.program
        setup = self.setup
        group_order = self.group_order

        if None not in witness:
            witness[None] = 0

        # Compute wire assignments for A, B, C, corresponding:
        # - A_values: witness[program.wires()[i].L]
        # - B_values: witness[program.wires()[i].R]
        # - C_values: witness[program.wires()[i].O]
        
        A_values = [Scalar(0) for _ in range(group_order)]
        B_values = [Scalar(0) for _ in range(group_order)]
        C_values = [Scalar(0) for _ in range(group_order)]
        
        for i in range(len(program.wires())):
            A_values[i] = (Scalar(witness[program.wires()[i].L]))
            B_values[i] = (Scalar(witness[program.wires()[i].R]))
            C_values[i] = (Scalar(witness[program.wires()[i].O]))
        
        # Construct A, B, C Lagrange interpolation polynomials for
        # A_values, B_values, C_values
        self.A = Polynomial(A_values, Basis.LAGRANGE)
        self.B = Polynomial(B_values, Basis.LAGRANGE)
        self.C = Polynomial(C_values, Basis.LAGRANGE)

        # Compute a_1, b_1, c_1 commitments to A, B, C polynomials
        a_1 = setup.commit(self.A)
        b_1 = setup.commit(self.B)
        c_1 = setup.commit(self.C)

        # Sanity check that witness fulfils gate constraints
        assert (
            self.A * self.pk.QL
            + self.B * self.pk.QR
            + self.A * self.B * self.pk.QM
            + self.C * self.pk.QO
            + self.PI
            + self.pk.QC
            == Polynomial([Scalar(0)] * group_order, Basis.LAGRANGE)
        )

        # Return a_1, b_1, c_1
        return Message1(a_1, b_1, c_1)

    def round_2(self) -> Message2:
        group_order = self.group_order
        setup = self.setup

        # Using A, B, C, values, and pk.S1, pk.S2, pk.S3, compute
        # Z_values for permutation grand product polynomial Z
        #
        # Note the convenience function:
        #       self.rlc(val1, val2) = val_1 + self.beta * val_2 + gamma
        
        # Z_values = [Scalar(1)]
        # roots_of_unity = Scalar.roots_of_unity(group_order)
        # for i in range(group_order):
        #     Z_values.append(
        #         Z_values[-1]
        #         * self.rlc(self.A.values[i], roots_of_unity[i])
        #         * self.rlc(self.B.values[i], 2 * roots_of_unity[i])
        #         * self.rlc(self.C.values[i], 3 * roots_of_unity[i])
        #         / self.rlc(self.A.values[i], self.pk.S1.values[i])
        #         / self.rlc(self.B.values[i], self.pk.S2.values[i])
        #         / self.rlc(self.C.values[i], self.pk.S3.values[i])
        #     )
        
        roots_of_unity = Scalar.roots_of_unity(group_order)
        
        Z_values = [Scalar(0) for _ in range(group_order + 1)]
        
        # the first value of the Z polynomial should be 1
        Z_values[0] = Scalar(1)
        
        for i in range(group_order):
            # F is RLC of a, b, c polynomial evaluated at original index and the original indices
            # indices are the ith root of unity
            # We make sure that domains do not overlap my multiplying by k_1 and k_2
            F_value = (
                self.rlc(self.A.values[i], roots_of_unity[i]) 
                * self.rlc(self.B.values[i], 2 * roots_of_unity[i])
                * self.rlc(self.C.values[i], 3 * roots_of_unity[i])
            )
            # G is RLC of a, b, c polynomials evaluated at the original index and the permutated index
            # This constitutes the "swap"
            G_value = (
                self.rlc(self.A.values[i], self.pk.S1.values[i])
                * self.rlc(self.B.values[i], self.pk.S2.values[i])
                * self.rlc(self.C.values[i], self.pk.S3.values[i])
            )
            # Z_value at i + 1 is the Z_value at i (previous Z value) multiplied with F/G
            Z_values[i + 1] = Z_values[i] * Scalar(F_value / G_value)

        # Check that the last term Z_n = 1
        assert Z_values.pop() == 1

        # Sanity-check that Z was computed correctly
        for i in range(group_order):
            assert (
                self.rlc(self.A.values[i], roots_of_unity[i])
                * self.rlc(self.B.values[i], 2 * roots_of_unity[i])
                * self.rlc(self.C.values[i], 3 * roots_of_unity[i])
            ) * Z_values[i] - (
                self.rlc(self.A.values[i], self.pk.S1.values[i])
                * self.rlc(self.B.values[i], self.pk.S2.values[i])
                * self.rlc(self.C.values[i], self.pk.S3.values[i])
            ) * Z_values[
                (i + 1) % group_order
            ] == 0

        # Construct Z, Lagrange interpolation polynomial for Z_values
        self.Z = Polynomial(Z_values, Basis.LAGRANGE)
        # Compute z_1 commitment to Z polynomial
        z_1 = setup.commit(self.Z)
        # Return z_1
        return Message2(z_1)

    def round_3(self) -> Message3:
        group_order = self.group_order
        setup = self.setup

        # Compute the quotient polynomial

        # List of roots of unity at 4x fineness, i.e. the powers of µ
        # where µ^(4n) = 1
        quarter_roots = Scalar.roots_of_unity(4 * group_order)

        # Using self.fft_expand, move A, B, C into coset extended Lagrange basis
        A_big = self.fft_expand(self.A)
        B_big = self.fft_expand(self.B)
        C_big = self.fft_expand(self.C)

        # Expand public inputs polynomial PI into coset extended Lagrange
        
        PI_big = self.fft_expand(self.PI)

        # Expand selector polynomials pk.QL, pk.QR, pk.QM, pk.QO, pk.QC
        # into the coset extended Lagrange basis
        
        QL_big = self.fft_expand(self.pk.QL)
        QR_big = self.fft_expand(self.pk.QR)
        QM_big = self.fft_expand(self.pk.QM)
        QO_big = self.fft_expand(self.pk.QO)
        QC_big = self.fft_expand(self.pk.QC)

        # Expand permutation grand product polynomial Z into coset extended
        # Lagrange basis
        
        Z_big = self.fft_expand(self.Z)

        # Expand shifted Z(ω) into coset extended Lagrange basis
        
        Z_shifted_big = Z_big.shift(4)

        # Expand permutation polynomials pk.S1, pk.S2, pk.S3 into coset
        # extended Lagrange basis

        S1_big = self.fft_expand(self.pk.S1)
        S2_big = self.fft_expand(self.pk.S2)
        S3_big = self.fft_expand(self.pk.S3)
        
        # Compute Z_H = X^N - 1, also in evaluation form in the coset
        
        # First we create a polynomial that evaluate to the roots of unity at the roots of unity
        # so it is like f(x) = x
        poly_eval_to_roots_of_unity = Polynomial(Scalar.roots_of_unity(group_order), basis=Basis.LAGRANGE)
        # Then expand into coset extended Lagrange basis
        nth = self.fft_expand(poly_eval_to_roots_of_unity)
        # Then we multiply this polynom with itself group_order - 1 times, hence getting X^N
        z_h = nth
        for _ in range(group_order - 1):
            z_h *= nth
        # then we substract 1 so z_h evaluates to 4 at the roots of unity
        z_h -= Scalar(1)
        self.z_h = z_h
        
        # Compute L0, the Lagrange basis polynomial that evaluates to 1 at x = 1 = ω^0
        # and 0 at other roots of unity
        # list_eval_to_1_at_w0 = [Scalar(1)] + [Scalar(0)] * group_order
        # self.L0 = Polynomial(list_eval_to_1_at_w0, basis=Basis.LAGRANGE)
        # L0_big = self.fft_expand(self.L0)
        
        self.L0 = Polynomial([Scalar(1)] + [Scalar(0)] * (group_order - 1), Basis.LAGRANGE)

        # Expand L0 into the coset extended Lagrange basis
        L0_big = self.fft_expand(
            Polynomial([Scalar(1)] + [Scalar(0)] * (group_order - 1), Basis.LAGRANGE)
        )

        # Compute the quotient polynomial (called T(x) in the paper)
        # It is only possible to construct this polynomial if the following
        # equations are true at all roots of unity {1, w ... w^(n-1)}:
        # 1. All gates are correct:
        #    A * QL + B * QR + A * B * QM + C * QO + PI + QC = 0
        #
        # 2. The permutation accumulator is valid:
        #    Z(wx) = Z(x) * (rlc of A, X, 1) * (rlc of B, 2X, 1) *
        #                   (rlc of C, 3X, 1) / (rlc of A, S1, 1) /
        #                   (rlc of B, S2, 1) / (rlc of C, S3, 1)
        #    rlc = random linear combination: term_1 + beta * term2 + gamma * term3
        #
        # 3. The permutation accumulator equals 1 at the start point
        #    (Z - 1) * L0 = 0
        #    L0 = Lagrange polynomial, equal at all roots of unity except 1
        
        correctGates = (
            A_big * B_big * QM_big 
            + A_big * QL_big 
            + B_big * QR_big 
            + C_big * QO_big 
            + PI_big + QC_big
            ) / z_h
        
        k_1 = Scalar(2)
        k_2 = Scalar(3)
        alpha = self.alpha
        fft_cofactor = self.fft_cofactor
        quarters = Polynomial(quarter_roots, Basis.LAGRANGE) * fft_cofactor
        
        # Since the expanded polynomial are in lagrange basis we can multiply it by the quarter roots of unity 
        # Otherwise we could do X = Polynomial([Scalar(0)] + [Scalar(1)] + [Scalar(0)] * (group_order - 2), Basis.MONOMIAL) 
        # and then use X_big = self.fft_expand(X) to create the X polynomial
        # or alternatively reuse poly_eval_to_roots_of_unity
        # the shift shows that z has been built up accumulateviley
        X = Polynomial([Scalar(0)] + [Scalar(1)] + [Scalar(0)] * (group_order - 2), Basis.MONOMIAL).fft()
        X_big = self.fft_expand(X)
        self.X_big = X_big
        
        permutationAcc = (
            self.rlc(A_big, X_big) * self.rlc(B_big, X_big * k_1) * self.rlc(C_big, X_big * k_2)
        ) * Z_big * alpha / z_h
        permutationAcc -= (
            self.rlc(A_big, S1_big) * self.rlc(B_big, S2_big) * self.rlc(C_big, S3_big)
        ) * Z_shifted_big * alpha / z_h
        
        result_of_acc = (Z_big - Scalar(1)) * L0_big * alpha * alpha / z_h
        
        QUOT_big = (
            correctGates
            + permutationAcc
            + result_of_acc
        )
        
        # print("Correct gates: ", self.expanded_evals_to_coeffs(correctGates).values[-group_order:])
        # print("permutationAcc: ", self.expanded_evals_to_coeffs(permutationAcc).values[-group_order:])
        # print("last term ", self.expanded_evals_to_coeffs(result_of_acc).values[-group_order:])
        # print("Left: ", self.expanded_evals_to_coeffs(QUOT_big).values[-group_order:])
        # print("Right: ", [0] * group_order)
        # print("Cofactor: ", self.fft_cofactor)
        # print("Roots of unity: ", Scalar.roots_of_unity(group_order))
        # print(" Roots of unity quarter: ", Scalar.roots_of_unity(group_order * 4))
        # print("Quarter * fft.cofactor: ", (quarters * fft_cofactor).values)
        # print("Group order: ", group_order)

        # Sanity check: QUOT has degree < 3n
        assert (
            self.expanded_evals_to_coeffs(QUOT_big).values[-group_order:]
            == [0] * group_order
        )
        print("Generated the quotient polynomial")
        
        QUOT_coeffs = QUOT_big.coset_extended_lagrange_to_coeffs(fft_cofactor)
        
        # print("QOUT_coeffs size: ", len(QUOT_coeffs.values))
        # print("QOUT_coeffs: ", QUOT_coeffs.values)

        # Split up T into T1, T2 and T3 (needed because T has degree 3n - 4, so is
        # too big for the trusted setup)
        # The last 8 coeffs of QOUT are 0s
        T1 = Polynomial(QUOT_coeffs.values[:group_order], Basis.MONOMIAL).fft()
        self.T1 = T1
        T2 = Polynomial(QUOT_coeffs.values[group_order:2*group_order], Basis.MONOMIAL).fft()
        self.T2 = T2
        T3 = Polynomial(QUOT_coeffs.values[2*group_order:3*group_order], Basis.MONOMIAL).fft()
        self.T3 = T3

        # Sanity check that we've computed T1, T2, T3 correctly
        assert (
            T1.barycentric_eval(fft_cofactor)
            + T2.barycentric_eval(fft_cofactor) * fft_cofactor**group_order
            + T3.barycentric_eval(fft_cofactor) * fft_cofactor ** (group_order * 2)
        ) == QUOT_big.values[0]

        print("Generated T1, T2, T3 polynomials")

        # Compute commitments t_lo_1, t_mid_1, t_hi_1 to T1, T2, T3 polynomials
        t_lo_1 = setup.commit(T1)
        t_mid_1 = setup.commit(T2)
        t_hi_1 = setup.commit(T3)

        # Return t_lo_1, t_mid_1, t_hi_1
        return Message3(t_lo_1, t_mid_1, t_hi_1)

    def round_4(self) -> Message4:
        # Compute evaluations to be used in constructing the linearization polynomial.
        
        # We only open these polynomials at zeta because this way 
        # we can calculate the final commitment for the linearization polynomial
        # These values (e.g. a-bar will act as replacement for the original polynomials)

        # Compute a_eval = A(zeta)
        a_eval = self.A.barycentric_eval(self.zeta)
        self.a_eval = a_eval 

        # Compute b_eval = B(zeta)
        b_eval = self.B.barycentric_eval(self.zeta)
        self.b_eval = b_eval 

        # Compute c_eval = C(zeta)
        c_eval = self.C.barycentric_eval(self.zeta)
        self.c_eval = c_eval 

        # Compute s1_eval = pk.S1(zeta)
        s1_eval = self.pk.S1.barycentric_eval(self.zeta)
        self.s1_eval = s1_eval
        # Compute s2_eval = pk.S2(zeta)
        s2_eval = self.pk.S2.barycentric_eval(self.zeta)
        self.s2_eval = s2_eval
        # Compute z_shifted_eval = Z(zeta * ω)
        z_shifted_eval = self.Z.barycentric_eval(self.zeta*Scalar.root_of_unity(self.group_order))
        self.z_shifted_eval = z_shifted_eval
        
        # Return a_eval, b_eval, c_eval, s1_eval, s2_eval, z_shifted_eval
        return Message4(a_eval, b_eval, c_eval, s1_eval, s2_eval, z_shifted_eval)

    def round_5(self) -> Message5:
        # Evaluate the Lagrange basis polynomial L0 at zeta
        L0_eval = self.L0.barycentric_eval(self.zeta)
        # Evaluate the vanishing polynomial Z_H(X) = X^n - 1 at zeta
        z_h_eval = self.zeta ** (self.group_order) - 1

        # Move T1, T2, T3 into the coset extended Lagrange basis
        T1_big = self.fft_expand(self.T1)
        T2_big = self.fft_expand(self.T2)
        T3_big = self.fft_expand(self.T3)
        # Move pk.QL, pk.QR, pk.QM, pk.QO, pk.QC into the coset extended Lagrange basis
        QL_big = self.fft_expand(self.pk.QL)
        QR_big = self.fft_expand(self.pk.QR)
        QM_big = self.fft_expand(self.pk.QM)
        QO_big = self.fft_expand(self.pk.QO)
        QC_big = self.fft_expand(self.pk.QC)
        # Move Z into the coset extended Lagrange basis
        Z_big = self.fft_expand(self.Z)
        # Move pk.S3 into the coset extended Lagrange basis
        S3_big = self.fft_expand(self.pk.S3)
        
        PI_eval = self.PI.barycentric_eval(self.zeta)

        # Compute the "linearization polynomial" R. This is a clever way to avoid
        # needing to provide evaluations of _all_ the polynomials that we are
        # checking an equation betweeen: instead, we can "skip" the first
        # multiplicand in each term. The idea is that we construct a
        # polynomial which is constructed to equal 0 at Z only if the equations
        # that we are checking are correct, and which the verifier can reconstruct
        # the KZG commitment to, and we provide proofs to verify that it actually
        # equals 0 at Z
        #
        # In order for the verifier to be able to reconstruct the commitment to R,
        # it has to be "linear" in the proof items, hence why we can only use each
        # proof item once; any further multiplicands in each term need to be
        # replaced with their evaluations at Z, which do still need to be provided
        
        R_big_gates_poly = (
            QM_big * self.a_eval * self.b_eval 
            + QL_big * self.a_eval 
            + QR_big * self.b_eval 
            + QO_big * self.c_eval 
            + PI_eval
            + QC_big
            )
        
        R_big_permutation = (
            (
                Z_big 
                * self.rlc(self.a_eval, self.zeta) 
                * self.rlc(self.b_eval, 2 * self.zeta) 
                * self.rlc(self.c_eval, 3 * self.zeta)
                - (
                    (S3_big * self.beta + self.c_eval + self.gamma) 
                    * self.rlc(self.a_eval, self.s1_eval) 
                    * self.rlc(self.b_eval, self.s2_eval) 
                    * self.z_shifted_eval
                )
            ) * self.alpha
        )
        
        R_big_permutation_acc = (
            ((Z_big - Scalar(1)) * L0_eval) * self.alpha ** 2
        )
        
        R_big_quo = (
            (
                T1_big 
                + T2_big * self.zeta ** self.group_order 
                + T3_big * self.zeta ** (self.group_order * 2)
            ) * z_h_eval
        )
        
        R_big = (
            R_big_gates_poly
            + R_big_permutation
            + R_big_permutation_acc
            - R_big_quo
        )
        # TODO: Refactor
        # R_big = (
        #     QM_big * self.a_eval * self.b_eval + QL_big * self.a_eval + QR_big * self.b_eval + QO_big * self.c_eval + self.PI.barycentric_eval(self.zeta) + QC_big
        #     + (
        #         Z_big * self.rlc(self.a_eval, self.zeta) * self.rlc(self.b_eval, 2 * self.zeta) * self.rlc(self.c_eval, 3 * self.zeta)
        #         - (S3_big * self.beta + self.c_eval + self.gamma) * self.rlc(self.a_eval, self.s1_eval) * self.rlc(self.b_eval, self.s2_eval) * self.z_shifted_eval
        #     ) * self.alpha
        #     + ((Z_big - Scalar(1)) * L0_eval) * self.alpha ** 2
        #     - (T1_big + T2_big * self.zeta ** self.group_order + T3_big * self.zeta ** (2 * self.group_order)) * z_h_eval
        # )
        
        # R_gates_poly_coeffs = self.expanded_evals_to_coeffs(R_big_gates_poly).values
        # R_gates_poly = Polynomial(R_gates_poly_coeffs[:self.group_order], Basis.MONOMIAL).fft()
        # print("R_gates: ", R_gates_poly.barycentric_eval(self.zeta))
        
        # R_permutation_coeffs = self.expanded_evals_to_coeffs(R_big_permutation).values
        # R_permutation_poly = Polynomial(R_permutation_coeffs[:self.group_order], Basis.MONOMIAL).fft()
        # print("R_permutation: ", R_permutation_poly.barycentric_eval(self.zeta))
        
        # R_permutation_acc_coeffs = self.expanded_evals_to_coeffs(R_big_permutation_acc).values
        # R_permutation_acc_poly = Polynomial(R_permutation_acc_coeffs[:self.group_order], Basis.MONOMIAL).fft()
        # print("R_permutation_acc: ", R_permutation_acc_poly.barycentric_eval(self.zeta))
        
        # R_quo_coeffs = self.expanded_evals_to_coeffs(R_big_quo).values
        # R_quo_poly = Polynomial(R_quo_coeffs[:self.group_order], Basis.MONOMIAL).fft()
        # print("R_quo: ", R_quo_poly.barycentric_eval(self.zeta))
        
        # print("z_h_eval: ", z_h_eval)
        # print("T1_big_eval: ", T1_big.barycentric_eval(self.zeta))
        # print("T2_big_eval: ", T2_big.barycentric_eval(self.zeta))
        # print("T3_big_eval: ", T3_big.barycentric_eval(self.zeta))
        
        R_coeffs = self.expanded_evals_to_coeffs(R_big).values
        # Degree test of R (max. 3n)
        assert R_coeffs[self.group_order:] == [0] * (self.group_order * 3)
        R = Polynomial(R_coeffs[:self.group_order], Basis.MONOMIAL).fft()
        # Commit to R
        
        com_R = self.setup.commit(R)

        # Sanity-check R
        assert R.barycentric_eval(self.zeta) == 0

        print("Generated linearization polynomial R")

        # Generate proof that W(z) = 0 and that the provided evaluations of
        # A, B, C, S1, S2 are correct

        # Move A, B, C into the coset extended Lagrange basis
        A_big = self.fft_expand(self.A)
        B_big = self.fft_expand(self.B)
        C_big = self.fft_expand(self.C)
        # Move pk.S1, pk.S2 into the coset extended Lagrange basis
        S1_big = self.fft_expand(self.pk.S1)
        S2_big = self.fft_expand(self.pk.S2)
        S3_big = self.fft_expand(self.pk.S3)
        
        # In the COSET EXTENDED LAGRANGE BASIS,
        # Construct W_Z = (
        #     R
        #   + v * (A - a_eval)
        #   + v**2 * (B - b_eval)
        #   + v**3 * (C - c_eval)
        #   + v**4 * (S1 - s1_eval)
        #   + v**5 * (S2 - s2_eval)
        # ) / (X - zeta)
        
        W_Z = (
            R_big
            + (A_big - self.a_eval) * self.v
            + (B_big - self.b_eval) * (self.v ** 2)
            + (C_big - self.c_eval) * (self.v ** 3)
            + (S1_big - self.s1_eval) * (self.v ** 4)
            + (S2_big - self.s2_eval) * (self.v ** 5)
        ) / (self.X_big - self.zeta)
        
        W_z_coeffs = self.expanded_evals_to_coeffs(W_Z).values

        # Check that degree of W_z is not greater than n
        assert W_z_coeffs[self.group_order:] == [0] * (self.group_order * 3)

        # Compute W_z_1 commitment to W_z
        
        W_z_1 = self.setup.commit(W_Z)

        # Generate proof that the provided evaluation of Z(z*w) is correct. This
        # awkwardly different term is needed because the permutation accumulator
        # polynomial Z is the one place where we have to check between adjacent
        # coordinates, and not just within one coordinate.
        # In other words: Compute W_zw = (Z - z_shifted_eval) / (X - zeta * ω)
        
        W_zw = (Z_big - self.z_shifted_eval) / (self.X_big - self.zeta * Scalar.root_of_unity(self.group_order))
        
        W_zw_coeffs = self.expanded_evals_to_coeffs(W_zw).values

        # Check that degree of W_z is not greater than n
        assert W_zw_coeffs[self.group_order:] == [0] * (self.group_order * 3)

        # Compute W_z_1 commitment to W_z
        
        W_zw_1 = self.setup.commit(W_zw)

        print("Generated final quotient witness polynomials")

        # Return W_z_1, W_zw_1
        return Message5(W_z_1, W_zw_1)

    def fft_expand(self, x: Polynomial):
        return x.to_coset_extended_lagrange(self.fft_cofactor)

    def expanded_evals_to_coeffs(self, x: Polynomial):
        return x.coset_extended_lagrange_to_coeffs(self.fft_cofactor)

    def rlc(self, term_1, term_2):
        return term_1 + term_2 * self.beta + self.gamma
