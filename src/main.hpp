/**

Note: Only 0-based index scheme is used in this program

**/

#ifndef __main_hpp__
#define __main_hpp__

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <cstring>
#include <vector>
#include <map>
#include <stdexcept>
#include <algorithm>

using std::vector;
using std::map;
using std::cout;
using std::endl;
using std::cin;
using std::runtime_error;
using std::max_element;
using std::exception;
using std::cerr;
using std::max;
using std::ofstream;
using std::swap;

#endif

/*
// Higher order element structure used in hierarchical FEM.
// 1-vertices, 2-edges, 3-quad faces, 4-tri faces, 5-interior

//         7------14-------6
//        /|              /|
//       / |             / |
//     15  |   25      13  |
//     /  19      22   /  18
//    /    |          /    |
//   4------12-------5     |
//   | 23  |   26    | 21  |
//   |     3------10-|-----2
//   |    /          |    /
//  16   /  20      17   /
//   | 11      24    |  9
//   | /             | /
//   |/              |/
//   0-------8-------1


//           5
//          /|\
//         / | \
//        /  |  \
//      11   |  10
//      /   14    \
//     /     |     \
//    /      |      \
//   3--------9------4
//   |  17   |  16   |
//   |       2       |
//   |      / \      |
//   |     /   \     |
//  12    / 15  \   13
//   |   8       7   |
//   |  /         \  |
//   | /           \ |
//   |/             \|
//   0-------6-------1

//  F3 (bottom tri face): 18, F4(top tri face): 19, I0: 20
//            3
//           /|\
//          / | \
//         /  |  \
//        9   |   8
//       /    |    \
//      /     |     \
//     /      7      \
//    2-------|5------1
//     \      |      /
//      \     |     /
//       \    |    /
//        6   |   4
//         \  |  /
//          \ | /
//           \|/
//            0
// Faces: F0:10, F1: 11, F2: 12, F3: 13, I0: 14

//      3-----6-----2
//      |           |
//      |           |
//      7     8     5
//      |           |
//      |           |
//      0-----4-----1


//      2
//      | \
//      |   \
//      5     4
//      |   6   \
//      |         \
//      0-----3----1


//
//	0-----2-----1
//
*/
/*
For hierarchical FEM, the term ``node" is used to as an abstraction for vertex, edge, face and interior of the element. (This is opposed to Lagrange basis)
The term ''mode'' denotes the basis functions of the FEM mesh. Also the interior of 2D elements are not considered as faces in our FEM structure.
*/
