from pathlib import Path

from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Dropout, Flatten
from tensorflow.keras.layers import Convolution3D, MaxPooling3D
from tensorflow.keras.callbacks import Callback

import tensorflow as tf




class AdaptiveLearningRate(Callback):

    def __init__(self, probe_interval=2000, accuracy_delta=0.001, reduction_factor=0.75, min_lr=1e-5):
        super().__init__()

        self.probe_interval = probe_interval
        self.accuracy_delta = accuracy_delta
        self.reduction_factor = reduction_factor
        self.min_lr = min_lr
        self.batch_num = 0
        self.last_accuracy = 0.0

    def on_train_batch_end(self, logs):

        if self.batch_num % self.probe_interval == 0:
            delta = logs["accuracy"] - self.last_accuracy
            self.last_accuracy = logs["accuracy"]
            self.model.reset_metrics()
            if delta < self.accuracy_delta:
                old_lr = self.model.optimizer.lr.read_value()
                new_lr = max(self.min_lr, old_lr * self.reduction_factor)
                self.model.optimizer.lr.assign(new_lr)

        self.batch_num += 1



def create_model(num_voxels=20):

    input_shape = (7, num_voxels, num_voxels, num_voxels)

    model = Sequential()
    model.add(
        Convolution3D(
            100,
            (3, 3, 3),
            activation="relu",
            input_shape=input_shape,
            data_format="channels_first",
        )
    )
    model.add(
        Convolution3D(200, (3, 3, 3), activation="relu", data_format="channels_first")
    )
    model.add(MaxPooling3D(pool_size=(2, 2, 2), data_format="channels_first"))
    model.add(
        Convolution3D(200, (2, 2, 2), activation="relu", data_format="channels_first")
    )
    model.add(
        Convolution3D(400, (2, 2, 2), activation="relu", data_format="channels_first")
    )
    model.add(MaxPooling3D(pool_size=(2, 2, 2), data_format="channels_first"))
    model.add(Flatten())
    model.add(Dense(1000, activation="relu"))
    model.add(Dropout(0.5))
    model.add(Dense(100, activation="relu"))
    model.add(Dropout(0.2))
    model.add(Dense(20, activation="softmax"))

    return model



def get_model(final_model_folder: Path, n_voxels, lr):

    if final_model_folder.exists():
        p_model = tf.keras.models.load_model(str(final_model_folder))
        p_model.optimizer.lr.assign(lr)
        print(f"--- loaded final model from {final_model_folder}", flush=True)
        return p_model

    p_model = create_model(n_voxels)

    optimizer = tf.keras.optimizers.SGD(learning_rate=lr, momentum=0.6)

    p_model.compile(
        loss="categorical_crossentropy", optimizer=optimizer, metrics=["accuracy"]
    )
    print("--- compiled model", flush=True)

    return p_model
