var CompressionPlugin = require('compression-webpack-plugin');
var CopyPlugin = require('copy-webpack-plugin');
var path = require('path');

var OUTPUT_DIR = './dist';

const fileName = 'chapter-map-component';
const libraryName = 'ChapterMapComponent';

function webpackConfig({ minimize }) {
  return {
    mode: process.env.WEBPACK_SERVE ? 'development' : 'production',
    entry: {
      [minimize ? `${fileName}.min` : fileName]: './src/index.js'
    },
    resolve: {
      extensions: ['.js', '.jsx']
    },
    output: {
      path: path.join(process.cwd(), OUTPUT_DIR),
      publicPath: '/dist',
      libraryTarget: 'umd',
      libraryExport: 'default',
      library: libraryName,
      filename: '[name].js',
      globalObject: '(typeof self !== "undefined" ? self : this)'
    },
    module: {
      rules: [
        {
          test: /\.jsx?$/,
          exclude: /node_modules/,
          use: { loader: 'babel-loader' }
        }
      ]
    },
    externals: {
      'prop-types': {
        root: 'PropTypes',
        commonjs: 'prop-types',
        commonjs2: 'prop-types',
        amd: 'prop-types'
      },
      react: {
        root: 'React',
        commonjs: 'react',
        commonjs2: 'react',
        amd: 'react'
      },
    },
    devtool: 'source-map',
    plugins: [
      new CompressionPlugin(),
      new CopyPlugin([
        { from: './src/*.json', to: '', flatten: true },
      ]),
    ],
    optimization: {
      noEmitOnErrors: true,
      minimize
    }
  };
}

module.exports = [
  webpackConfig({ minimize: true}),
  webpackConfig({ minimize: false })
];
