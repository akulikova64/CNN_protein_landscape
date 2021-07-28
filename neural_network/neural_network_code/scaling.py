from pathlib import Path

import numpy as np




class UpdatingScaler:
    """A class capable of scaling a list of boxes
    and updating the average and standard deviation
    with the values it's scaling
    """

    def __init__(self, m, m2, count):
        """Constructor.
        :param m: the current average value
        :type m: tensor of float
        :param m2: a value related to the standard deviation
        :type m2: tensor of float
        :param count: the number of arrays used so far
        :type count: int
        The input arguments for this class are the output of function
        :meth:`online_generator`. They are the values described by
        `Welford's online algorithm <https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance>`_
        for computing mean and standard deviation.
        """
        self.m = m
        self.m2 = m2
        self.count = count


    def scale(self, new_vals, update=False):
        """Scale the boxes and potentially update the
        mean and standard deviation
        :param new_vals: the boxes to scale
        :type new_vals: a 5d tensor
        :param update: a flag denoting that the mean and standard deviation should be updated, defaults to False
        :type update: bool, optional
        """

        if update:
            self.m, self.m2, self.count = online_arraylike(
                new_vals, self.m, self.m2, self.count
            )

        ave, std = self.get_ave_std()
        std[std == 0] = np.float64(1)

        for i, v in enumerate(new_vals):
            new_vals[i] = (v - ave) / std


    def get_ave_std(self):
        return self.m, np.sqrt(self.m2 / float(self.count))



def _update(prev_mean, prev_m2, new_val, count):
    count += 1
    delta = new_val - prev_mean
    new_mean = prev_mean + delta / float(count)
    delta2 = new_val - new_mean
    m2 = prev_m2 + delta * delta2
    return count, new_mean, m2



def online_arraylike(vals, m, m2, count):
    for v in vals:
        count, m, m2 = _update(m, m2, v, count)
    return m, m2, count



def load_scaler(output_dir: Path):
    scaler_folder = output_dir / "scaler"
    if not scaler_folder.exists():
        scaler_folder.mkdir()
        return None
    else:
        cnt_file = scaler_folder / "count.txt"
        if not cnt_file.exists():
            return None
        with cnt_file.open("r") as f:
            count = f.read()
        m = np.load(scaler_folder / "mean.npy")
        m2 = np.load(scaler_folder / "m2.npy")

        print(f"--- loaded previous scaler")

        return UpdatingScaler(m, m2, count)



def save_scaler(scaler: UpdatingScaler, output_dir: Path):
    scaler_folder = output_dir / "scaler"
    if not scaler_folder.exists():
        scaler_folder.mkdir()

    np.save(scaler_folder / "mean", scaler.m)
    np.save(scaler_folder / "m2", scaler.m2)
    cnt_file = scaler_folder / "count.txt"
    with cnt_file.open("w") as f:
        f.write(str(scaler.count))
        