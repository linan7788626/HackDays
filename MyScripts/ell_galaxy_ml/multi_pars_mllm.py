import matplotlib.pyplot as plt
import numpy as np
import astropy.io.fits as pyfits
import subprocess as sp
import sys
import os
import cPickle

import theano
import theano.tensor as T
import lasagne
from lasagne.layers import Conv2DLayer as conv
from lasagne.layers import MaxPool2DLayer as pool
from lasagne.layers import DenseLayer as dense
from lasagne.nonlinearities import rectify
from lasagne.regularization import l2, regularize_layer_params


def ReadModelFromFile(model, filename):
    if not os.path.isfile(filename):
        print (filename, "not exist")
        sys.exit(-1)
    f = open(filename, 'r')
    data = cPickle.load(f)
    lasagne.layers.set_all_param_values(model, data)


def WriteModelToFile(model, filename):
    data = lasagne.layers.get_all_param_values(model)
    with open(filename, 'w') as f:
        cPickle.dump(data, f)


def BuildModel(x, img_size=500, num_classes=3):
    filter_size = [3,3,3,3,3,3]      # 3x3 or 5x5
    pool_size =   [1,1,2,2,2,2]      # pooling

    num_filters = [16, 16, 32, 32, 32, 32]

    n_layers = len(filter_size)
    if any( [len(pool_size) != n_layers, len(num_filters)!=n_layers] ):
        raise ValueError('num of layers should be consistent')

    n_hiddens = [10, 10]      #[50, 50]        public
    n_hiddens_self = [10, 10] #[100, 100]      private

    model = lasagne.layers.InputLayer(shape=(None, 1, img_size, img_size), input_var=x)

    for i in xrange(n_layers):
        model = conv( model, num_filters=num_filters[i], filter_size=(filter_size[i], filter_size[i]),
                     nonlinearity=rectify, W=lasagne.init.GlorotUniform(), pad='same')
        if pool_size[i] > 1:
            model = pool(model, pool_size=(pool_size[i],pool_size[i]))

    for i in xrange( len(n_hiddens) ):
        model = dense(model, num_units = n_hiddens[i], nonlinearity=rectify)

    l_out = []   #FINAL MODEL

    for i in xrange(num_classes):
        tmp = model
        for j in xrange( len(n_hiddens_self) ):
            tmp = dense(tmp, num_units=n_hiddens_self[j], nonlinearity=rectify)
        tmp = dense(tmp, num_units=1, nonlinearity=rectify)
        l_out.append(tmp)
        tmp = None
    return l_out


def Fit(train_model, X_train, Y_train, N_EPOCH, batch_train,
        valid_model, X_valid, Y_valid, batch_valid, write_to='./'):

    if not os.path.isdir:
        raise ValueError('path not exists: ' + write_to)
    n_round = X_train.shape[0]/ batch_train + (X_train.shape[0] % batch_train > 0)

    loss_train = []
    loss_valid = []

    for e in xrange(N_EPOCH):
        sys.stdout.write("\r%d / %d" % (e+1, N_EPOCH ) )
        sys.stdout.flush()

        loss_epoch = []
        for idx in xrange( n_round ):
            start_index = idx * batch_train
            end_index = (idx+1) * batch_train
            if end_index > X_train.shape[0]-1:
                end_index = X_train.shape[0]-1
            if start_index >= end_index:
                continue
            tmp = train_model(X_train[start_index:end_index], Y_train[start_index:end_index,0],
                               Y_train[start_index:end_index, 1], Y_train[start_index:end_index, 2])
            loss_epoch.append( np.array( tmp ) )
        loss_train.append( np.mean( np.row_stack(loss_epoch), 0 ) )

        #WriteModelToFile(l_out, os.path.join(write_to, str(e+1)+'.pkl') )

        n_round_valid = X_valid.shape[0]/ batch_valid + (X_valid.shape[0] % batch_valid > 0)

        loss_epoch = []
        for idx in xrange( n_round_valid ):
            start_index = idx * batch_valid
            end_index = (idx+1) * batch_valid
            if end_index > X_valid.shape[0]-1:
                end_index = X_valid.shape[0]-1
            if start_index >= end_index:
                continue
            tmp = valid_model(X_valid[start_index:end_index], Y_valid[start_index:end_index,0],
                               Y_valid[start_index:end_index, 1], Y_valid[start_index:end_index, 2])
            loss_epoch.append( np.array( tmp ) )
        loss_valid.append( np.mean( np.row_stack(loss_epoch), 0 ) )

    loss_train = np.array(loss_train)
    loss_valid = np.array(loss_valid)

    return loss_train, loss_valid


def main():
    fits_dir = "./fits_outputs_0/"
    n_images = 100
    frat = 0.8

    cmd = "ls " + fits_dir
    fits_list = sp.check_output(cmd,shell=True)
    fits_files = fits_list.split('\n')[:-1]

    parameters = []
    images = []
    for i in fits_files[:n_images]:
        ftmp = i.split('_')
        parameters.append(np.array(ftmp[1:-1]).astype("float32"))
        fname = fits_dir + i
        itmp = pyfits.getdata(fname)
        images.append(np.array(itmp).astype("float32"))

    imgs = np.array(images)
    pars = np.array(parameters)
    X = imgs.reshape(imgs.shape[0], 1, imgs.shape[1], imgs.shape[2])

    Y = pars[:, [2, 3, 4]]

    Y[:,0] = (Y[:,0] - 0.5) / 0.5 * 100.0
    Y[:,1] = (Y[:,1] - 1.0) / 3.0 * 100.0
    Y[:,2] = (Y[:,2]) / 360.0 * 100.0

    train_index = range(0, int(frat * n_images))
    valid_index = range(int(frat * n_images), n_images)

    train_x = X[train_index]
    train_y = Y[train_index]
    valid_x = X[valid_index]
    valid_y = Y[valid_index]
#----------------------------------------------------------------------------
    num_classes = 3
    img_size = 500
    lambda_reg = 0.0001

    x = T.tensor4('input image')
    y = []
    for i in xrange(num_classes):
        y.append( T.vector() )

    l_out = BuildModel(x, img_size, num_classes)   #model name

    pred_train = lasagne.layers.get_output(l_out, deterministic=False)
    pred_valid = lasagne.layers.get_output(l_out, deterministic=True)

    params = lasagne.layers.get_all_params(l_out, trainable=True, regularizable=True)

    mse_loss_train = []
    for i in xrange(num_classes):
        mse_loss_train.append(lasagne.objectives.squared_error(y[i], pred_train[i]).mean())

    reg_loss = regularize_layer_params(l_out, l2) * lambda_reg
    loss = T.mean(mse_loss_train) + reg_loss
    updates = lasagne.updates.adam(loss, params)
    #train_model = theano.function([x] + y, [loss] + mse_loss_train + [reg_loss], updates=updates)
    train_model = theano.function([x] + y, [loss] + mse_loss_train, updates=updates)

    mse_loss_valid = []
    for i in xrange(num_classes):
        mse_loss_valid.append( lasagne.objectives.squared_error(y[i], pred_valid[i]).mean() )
    valid_model = theano.function([x] + y, [loss] + mse_loss_train)
#----------------------------------------------------------------------------
    N_EPOCH = 100
    batch_train = 1
    batch_valid = 1
    models_dir = "./models_1000_0/"

    loss_train, loss_valid = Fit(train_model, train_x, train_y, N_EPOCH, batch_train,
                                 valid_model, valid_x, valid_y, batch_valid,
                                 write_to=models_dir)

    #loss_train.tofile(models_dir+"loss_train.bin", format="float32")
    #loss_valid.tofile(models_dir+"loss_valid.bin", format="float32")
#----------------------------------------------------------------------------
    fig = plt.figure( figsize=(10, 4) )
    fig.add_subplot(121)
    plt.plot(loss_train, 'o--')
    plt.title('Training')
    plt.xlabel("Epoch")
    plt.ylabel("Loss")
    plt.legend(("Total", "Ellipticity", "Einstein R", "Orientation"),
               loc='upper right', shadow=True, fontsize=14)

    fig.add_subplot(122)
    plt.plot(loss_valid, 'o--')
    plt.title('Valid')
    plt.xlabel("Epoch")
    plt.ylabel("Loss")
    plt.legend(("Total", "Ellipticity", "Einstein R", "Orientation"),
               loc='upper left', shadow=True, fontsize=14)
#----------------------------------------------------------------------------
    f_pred = theano.function([x], pred_train)
    y_pred = []

    y_pred = f_pred(X[:])
    y_pred = np.array(y_pred)[:,:,0]
#----------------------------------------------------------------------------
    fig = plt.figure( figsize= (15, 4) )

    fig.add_subplot(131)
    plt.title('Ellipticity')
    plt.plot(y_pred[0, train_index], Y[train_index,0], 'o')
    plt.plot(np.linspace(0,149,150), np.linspace(0,149,150), "r-")
    plt.axis([0, 150, 0, 150])

    fig.add_subplot(132)
    plt.title('Einstein R')
    plt.plot(y_pred[1, train_index], Y[train_index,1], 'o')
    plt.plot(np.linspace(0,149,150), np.linspace(0,149,150), "r-")
    plt.axis([0, 150, 0, 150])

    fig.add_subplot(133)
    plt.title('Orientations')
    plt.plot(y_pred[2, train_index], Y[train_index,2], 'o')
    plt.plot(np.linspace(0,149,150), np.linspace(0,149,150), "r-")
    plt.axis([0, 150, 0, 150])
#----------------------------------------------------------------------------
    fig = plt.figure( figsize= (15, 4) )

    fig.add_subplot(131)
    plt.plot(y_pred[0, valid_index], Y[valid_index,0], 'o')
    plt.plot(np.linspace(0,149,150), np.linspace(0,149,150), "r-")
    plt.axis([0, 150, 0, 150])

    fig.add_subplot(132)
    plt.plot(y_pred[1, valid_index], Y[valid_index,1], 'o')
    plt.plot(np.linspace(0,149,150), np.linspace(0,149,150), "r-")
    plt.axis([0, 150, 0, 150])

    fig.add_subplot(133)
    plt.plot(y_pred[2, valid_index], Y[valid_index,2], 'o')
    plt.plot(np.linspace(0,149,150), np.linspace(0,149,150), "r-")
    plt.axis([0, 150, 0, 150])
#----------------------------------------------------------------------------
    return 0

if __name__ == '__main__':
    main()
    plt.show()
