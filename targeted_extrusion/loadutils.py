import numpy as np
import pandas as pd
import h5py
from copy import deepcopy


class leg(object):
    def __init__(
        self, pos, attrs={"stalled": False, "CTCF": False, "in solution": False}
    ):
        """
        A leg has two important attribues: pos (positions) and attrs (a custom list of attributes)
        """
        self.pos = pos
        self.attrs = dict(attrs)


class cohesin(object):
    """
    A cohesin class provides fast access to attributes and positions


    cohesin.left is a left leg of cohesin, cohesin.right is a right leg
    cohesin[-1] is also a left leg and cohesin[1] is a right leg

    Also, cohesin.any("myattr") is True if myattr==True in at least one leg
    cohesin.all("myattr") is if myattr=True in both legs
    """

    def __init__(self, leg1, leg2):
        self.loading_point = deepcopy(leg1)
        self.left = leg1
        self.right = leg2

    def any(self, attr):
        return self.left.attrs[attr] or self.right.attrs[attr]

    def all(self, attr):
        return self.left.attrs[attr] and self.right.attrs[attr]

    def __getitem__(self, item):
        if item == -1:
            return self.left
        elif item == 1:
            return self.right
        elif item == 0:
            return self.loading_point
        else:
            raise ValueError()


def unloadProb(cohesin, args):
    """
    Defines unload probability based on a state of cohesin
    """
    if cohesin.any("stalled"):
        # if one side is stalled, we have different unloading probability
        # Note that here we define stalled cohesins as those stalled not at CTCFs
        return 1 / args["LIFETIME_STALLED"]
    # otherwise we are just simply unloading
    return 1 / args["LIFETIME"]


def loadOne(cohesins, occupied, args):
    """
    A function to load one cohesin
    """
    while True:
        a = np.random.randint(args["N"])
        if (occupied[a] == 0) and (occupied[a + 1] == 0):
            occupied[a] = 1
            occupied[a + 1] = 1
            cohesin(leg(a), leg(a + 1))
            for side in [-1, 1]:
                c[side].attrs["loaded"] = True
            cohesins.append(c)
            break


def loadOneWithProbs(cohesins, occupied, args, position=None):
    """
    A function to load one cohesin given the track of probabilities
    """

    # If nothing is loadable, let cohesins stay in the solution:
    if np.all(np.array(occupied) == 1) or np.all(
        np.array(args["probsLoading"]) * (1 - np.array(occupied)) == 0
    ):
        c = cohesin(leg(np.nan), leg(np.nan))
        for side in [-1, 1]:
            c[side].attrs["in solution"] = True
            c[side].attrs["stalled"] = False
            c[side].attrs["loaded"] = False
        if position is None:
            cohesins.append(c)
        else:
            cohesins[position] = c
        # raise ValueError('Unable to load!')

    else:
        ntrials = 0
        while ntrials < 3:
            a = np.random.choice(np.arange(args["N"]), size=1, p=args["probsLoading"])[
                0
            ]
            if (occupied[a] == 0) and (occupied[a + 1] == 0):
                occupied[a] = 1
                occupied[a + 1] = 1
                c = cohesin(leg(a), leg(a + 1))
                for side in [-1, 1]:
                    c[side].attrs["loaded"] = True
                    c[side].attrs["stalled"] = False
                if position is None:
                    cohesins.append(c)
                else:
                    cohesins[position] = c
                break
            ntrials += 1

        # If nothing is loadable, let cohesins stay in the solution:
        else:
            c = cohesin(leg(-1), leg(-1))
            for side in [-1, 1]:
                c[side].attrs["in solution"] = True
                c[side].attrs["loaded"] = False
                c[side].attrs["stalled"] = False
            if position is None:
                cohesins.append(c)
            else:
                cohesins[position] = c


def capture(cohesin, occupied, args):
    """
    We are describing CTCF capture here.
    This function is specific to this particular project, and
    users are encouraged to write functions like this

    Note the for-loop over left/right sites below, and using cohesin[side]
    to get left/right leg.

    Also note how I made ctcfCapture a dict with -1 coding for left side, and 1 for right side
    and ctcfCapture are dicts as well: keys are locations, and values are probabilities of capture
    """
    for side in [1, -1]:
        # get probability of capture or otherwise it is 0
        if np.random.random() < args["ctcfCapture"][side].get(cohesin[side].pos, 0):
            cohesin[side].attrs["CTCF"] = True  # captured a cohesin at CTCF
    return cohesin


def release(cohesin, occupied, args):
    """
    AN opposite to capture - releasing cohesins from CTCF
    """

    if not cohesin.any("CTCF"):
        return cohesin  # no CTCF: no release necessary

    # attempting to release either side
    for side in [-1, 1]:
        if (
            np.random.random() < args["ctcfRelease"][side].get(cohesin[side].pos, 0)
        ) and (cohesin[side].attrs["CTCF"]):
            cohesin[side].attrs["CTCF"] = False
    return cohesin


def translocate(cohesins, occupied, args):
    """
    This function describes everything that happens with cohesins -
    loading/unloading them and stalling against each other

    It relies on the functions defined above: unload probability, capture/release.
    """

    # Change the status of all loaded cohesins to moving:
    for i in range(len(cohesins)):
        cohesin = cohesins[i]
        if cohesin.any("loaded"):
            for leg in [-1, 1]:
                cohesin[leg].attrs["loaded"] = False

    # first we try to unload cohesins and free the matching occupied sites
    i = 0
    last_tocheck = len(cohesins)
    while i < last_tocheck:
        if not cohesins[i].any("in solution"):
            prob = unloadProb(cohesins[i], args)
            if np.random.random() < prob:
                occupied[cohesins[i].left.pos] = 0
                occupied[cohesins[i].right.pos] = 0
                # del cohesins[i]
                loadOneWithProbs(
                    cohesins, occupied, args, i
                )  # Immediately load cohesin when it detaches
                # last_tocheck += -1
            # else:
            # i += 1
        # else:
        i += 1

    # # then we try to capture and release them by CTCF sites
    for i in range(len(cohesins)):
        #         print(f'capture-check {i}')
        cohesins[i] = capture(cohesins[i], occupied, args)
        cohesins[i] = release(cohesins[i], occupied, args)

    # Try to load cohesins from the solution
    i = 0
    last_tocheck = len(cohesins)
    while i < last_tocheck:
        if cohesins[i].any("in solution"):
            # del cohesins[i]
            loadOneWithProbs(cohesins, occupied, args, i)
            # last_tocheck += -1
        # else:
        i += 1

    # finally we translocate, and mark stalled cohesins because the unloadProb needs this
    for i in range(len(cohesins)):
        cohesin = cohesins[i]
        if cohesin.any("in solution"):
            continue
        if cohesin.any("loaded"):
            continue
        else:
            for leg in [-1, 1]:
                try:
                    if not cohesin[leg].attrs["CTCF"]:
                        # cohesins that are not at CTCFs and cannot move are labeled as stalled
                        if occupied[cohesin[leg].pos + leg] != 0:
                            cohesin[leg].attrs["stalled"] = True
                        else:
                            cohesin[leg].attrs["stalled"] = False
                            occupied[cohesin[leg].pos] = 0
                            occupied[cohesin[leg].pos + leg] = 1
                            cohesin[leg].pos += leg
                except Exception as e:
                    print(cohesin[leg].attrs, cohesin[leg].pos, e)
            cohesins[i] = cohesin


def color(cohesins, args):
    "A helper function that converts a list of cohesins to an array colored by cohesin state"

    def state(attrs):
        if attrs["stalled"]:
            return 2
        if attrs["CTCF"]:
            return 3
        return 1

    ar = np.zeros(args["N"])
    for i in cohesins:
        ar[i.left.pos] = state(i.left.attrs)
        ar[i.right.pos] = state(i.right.attrs)
    return ar


def color_targeted(cohesins, args):
    "A helper function that converts a list of cohesins to an array colored by cohesin state"

    def state(attrs):
        if attrs["stalled"]:
            return 5
        if attrs["CTCF"]:
            return 6
        return 4

    ar = np.zeros(args["N"])
    for i in cohesins:
        ar[i.left.pos] = state(i.left.attrs)
        ar[i.right.pos] = state(i.right.attrs)
    return ar


# Additional functions:
def SparseToMatrix(sparse_contacts, min_coord, max_coord):
    n = max_coord - min_coord
    matrix = np.zeros([n, n])
    for left, right in sparse_contacts:
        if left >= min_coord and right < max_coord:
            matrix[left - min_coord, right - min_coord] += 1

    return matrix


def makeMatrix(cohesins, args):
    interaction_matrix = np.zeros([args["N"], args["N"]])

    for i in cohesins:
        interaction_matrix[i.left.pos, i.right.pos] += 1  # state(i.left.attrs)

    return interaction_matrix


def makeSparseContacts(cohesins, args):
    contacts = []
    for i in cohesins:
        contacts.append([i.left.pos, i.right.pos])
    return contacts
