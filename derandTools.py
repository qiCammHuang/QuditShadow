import math


def derandomized_classical_shadow_qudit(all_observables, measurement_budget, system_size, *, dim=3, eta_default=0.9,
                                        weight=None, local_bias=None):
    """
        Implementation of derandomized classical shadow for qudits (qutrits as default)
        A generalized version of codes originally by Hsin-Yuan Huang (https://momohuang.github.io/)

        all_observables: a list of orthonormal Hermitian observables which, in our case, is taken as a set of
        generalized Gell-Mann Matrices indexed by 0, 1, ..., (dim^2-1). And is actually stored as a list of tuples of
        the form (1, position) or (4, position) or (7, position)

        measurement_budget: int for the amount of measurement budget

        system_size: int for the number of qudits in the system

        dim: int for qudit-dimensions. 2 for qubits, 3 for qutrits

        weight: None or a list of coefficients for each observable
                None -- neglect this parameter
                a list -- modify the number of measurements for each observable by the corresponding weight

        local_bias: None or a list of lists indicating an a priori probability distribution bias as opposed to uniform
                None -- automatically generate uniform sampling processes
                a list of lists -- modify the a priori prob. dist. for each qudits, with the encompassed lists
                                   indicating actionable advice from LBCS optimization
                    local_bias[qudit index][ggm index]
    """

    if weight is None:
        weight = [1.0] * len(all_observables)
    assert (len(weight)) == len(all_observables)

    if local_bias is None:
        local_bias = [[1.0] * (dim ** 2) for _ in range(system_size)]
    assert (len(local_bias) == system_size)

    for qudit in range(system_size):
        local_bias[qudit][0] = 0.0

    local_bias_normalizer = []
    for qudit_i, prob in enumerate(local_bias):
        assert (len(prob) == dim ** 2)
        local_bias_normalizer.append(sum(prob) - prob[0])

    sum_log_value = 0
    sum_cnt = 0

    def success_average(num_of_measurements_so_far):
        # an indicator of average hit level per observable, biased upon weight.
        result = 0
        for i in range(len(all_observables)):
            result += num_of_measurements_so_far[i] / weight[i]
        result /= len(all_observables)
        return result

    def cost_function_qudit(num_of_measurements_so_far, num_of_matches_needed_in_this_round, derand_progress):
        # derand_progress: an integer, number of fully derandomized basis samples, correspondingly (m-1) in paper
        #                  0 for the first round, (M-1) for the last round within measurement budget.

        eta = eta_default  # hyper-param subject to change
        nu = 1 - math.exp(- eta / 2)

        nonlocal sum_log_value
        nonlocal sum_cnt

        cost = 0

        for i, zipitem in enumerate(zip(num_of_measurements_so_far, num_of_matches_needed_in_this_round)):
            current_observable = all_observables[i]
            measurement_so_far, matches_needed = zipitem

            # I altered the idea of "success" here - qi huang
            success_standard = measurement_budget
            if derand_progress >= 0.1 * measurement_budget:
                success_standard = success_average(num_of_measurements_so_far)

            # introduced to rescale cost function, alleviate the 'stuck' issue for large M's
            shift = eta * success_average(num_of_measurements_so_far) / 2

            if num_of_measurements_so_far[i] >= math.floor(weight[i] * success_standard) + 1:
                # successful observables are not considered in evaluting the cost function
                continue
            if system_size < matches_needed:
                V = eta / 2 * measurement_so_far
            else:
                biased_prob = 1.0

                cnt_flag = 0
                for ggm_unit in reversed(current_observable):
                    if cnt_flag > matches_needed - 0.5:
                        break
                    else:
                        biased_prob *= local_bias[ggm_unit[1]][ggm_unit[0]] / local_bias_normalizer[ggm_unit[1]]
                        cnt_flag += 1

                V = eta / 2 * measurement_so_far - math.log(1 - nu * biased_prob)
            cost += math.exp(-V / weight[i] + shift)

            sum_log_value += V / weight[i]
            sum_cnt += 1

        return cost

    def match_up_qudit(qudit_i, dice_roll_ggm, single_observable):
        for ggm, pos in single_observable:
            if pos != qudit_i:
                continue
            else:
                if ggm != dice_roll_ggm:
                    return -1
                else:
                    return 1
        return 0

    num_of_measurements_so_far = [0] * len(all_observables)
    measurement_procedure = []

    # debug 'stuck' issue
    # stuck_flag = 0

    for repetition in range(measurement_budget):
        # a single round of parallel measurement over "system_size" number of qudits
        num_of_matches_needed_in_this_round = [len(P) for P in all_observables]  # Pauli-weights (# of non-identities)
        single_round_measurement = []

        shift = sum_log_value / sum_cnt if sum_cnt > 0 else 0;
        sum_log_value = 0.0
        sum_cnt = 0

        for qudit_i in range(system_size):
            cost_of_outcomes = []
            list_of_ggm = []
            for ggm in range(1, dim ** 2):
                cost_of_outcomes.append((str(ggm), 0.0))
                list_of_ggm.append(ggm)
            cost_of_outcomes = dict(cost_of_outcomes)

            for dice_roll_ggm in list_of_ggm:
                # Roll the dice for GGM (as opposed of Pauli matrices)
                for i, single_observable in enumerate(all_observables):
                    result = match_up_qudit(qudit_i, dice_roll_ggm, single_observable)
                    if result == -1:
                        num_of_matches_needed_in_this_round[i] += 100 * (system_size + 10)  # impossible to measure
                    if result == 1:
                        num_of_matches_needed_in_this_round[i] -= 1  # match up one Generalized Gell-Mann Matrix

                # Record cost function outputs
                cost_of_outcomes[str(dice_roll_ggm)] = cost_function_qudit(num_of_measurements_so_far,
                                                                           num_of_matches_needed_in_this_round,
                                                                           derand_progress=repetition)

                '''
                # debug 'stuck' issue
                if stuck_flag > 0.5:
                    print(cost_of_outcomes)
                    print(success_average(num_of_measurements_so_far))
                    for i in range(len(all_observables)):
                        print(num_of_measurements_so_far[i] / weight[i])
                '''

                # Revert the dice roll
                for i, single_observable in enumerate(all_observables):
                    result = match_up_qudit(qudit_i, dice_roll_ggm, single_observable)
                    if result == -1:
                        num_of_matches_needed_in_this_round[i] -= 100 * (system_size + 10)
                    if result == 1:
                        num_of_matches_needed_in_this_round[i] += 1

            for dice_roll_ggm in list_of_ggm:
                if min(cost_of_outcomes.values()) < cost_of_outcomes[str(dice_roll_ggm)]:
                    continue
                # The best dice roll outcome will come to this line
                single_round_measurement.append(dice_roll_ggm)
                for i, single_observable in enumerate(all_observables):
                    result = match_up_qudit(qudit_i, dice_roll_ggm, single_observable)
                    if result == -1:
                        num_of_matches_needed_in_this_round[i] += 100 * (system_size + 10)  # impossible to measure
                    if result == 1:
                        num_of_matches_needed_in_this_round[i] -= 1  # match up one Generalized Gell-Mann Matrix
                break

        measurement_procedure.append(single_round_measurement)

        '''
        # introduced to debug the 'stuck' issue
        if single_round_measurement == [1] * system_size:
            stuck_flag = 1
        else:
            stuck_flag = 0
        '''

        for i, single_observable in enumerate(all_observables):
            if num_of_matches_needed_in_this_round[i] == 0:  # finished measuring all qudits
                num_of_measurements_so_far[i] += 1

    return measurement_procedure
