# This file is part of PG.
#
# PG is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# PG is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with PG.  If not, see <http://www.gnu.org/licenses/>.

from pybindgen import *
import sys


def import_rbd_types(mod):
  mod.add_class('MultiBody', foreign_cpp_namespace='rbd', import_from_module='rbdyn')


def import_sva_types(mod):
  mod.add_class('MotionVecd', foreign_cpp_namespace='sva', import_from_module='spacevecalg')
  mod.add_class('ForceVecd', foreign_cpp_namespace='sva', import_from_module='spacevecalg')
  mod.add_class('RBInertiad', foreign_cpp_namespace='sva', import_from_module='spacevecalg')
  mod.add_class('ABInertiad', foreign_cpp_namespace='sva', import_from_module='spacevecalg')
  mod.add_class('PTransformd', foreign_cpp_namespace='sva', import_from_module='spacevecalg')



def import_eigen3_types(mod):
  mod.add_class('Vector3d', foreign_cpp_namespace='Eigen', import_from_module='eigen3')
  mod.add_class('Vector6d', foreign_cpp_namespace='Eigen', import_from_module='eigen3')

  mod.add_class('Matrix3d', foreign_cpp_namespace='Eigen', import_from_module='eigen3')
  mod.add_class('Matrix6d', foreign_cpp_namespace='Eigen', import_from_module='eigen3')

  mod.add_class('MatrixXd', foreign_cpp_namespace='Eigen', import_from_module='eigen3')
  mod.add_class('VectorXd', foreign_cpp_namespace='Eigen', import_from_module='eigen3')

  mod.add_class('Quaterniond', foreign_cpp_namespace='Eigen', import_from_module='eigen3')



def build_pg(pg):
  pgSolver = pg.add_class('PostureGenerator', template_parameters=['pg::eigen_ad'], custom_name='PostureGenerator')

  fixedContact = pg.add_struct('FixedContact')

  # build list type
  pg.add_container('std::vector<pg::FixedContact>',
                      'pg::FixedContact', 'vector')

  # PostureGenerator
  pgSolver.add_constructor([param('const rbd::MultiBody&', 'mb')])

  pgSolver.add_method('fixedContacts', None, [param('std::vector<pg::FixedContact>', 'contacts')])
  pgSolver.add_method('run', retval('bool'), [])


  # FixedContact
  fixedContact.add_constructor([])

  fixedContact.add_instance_attribute('bodyId', 'int')
  fixedContact.add_instance_attribute('pos', 'Eigen::Vector3d')



if __name__ == '__main__':
  if len(sys.argv) < 2:
    sys.exit(1)

  pg = Module('_pg', cpp_namespace='::pg')
  pg.add_include('<EigenAutoDiffScalar.h>')
  pg.add_include('<PostureGenerator.h>')

  dom_ex = pg.add_exception('std::domain_error', foreign_cpp_namespace=' ',
                               message_rvalue='%(EXC)s.what()')
  out_ex = pg.add_exception('std::out_of_range', foreign_cpp_namespace=' ',
                               message_rvalue='%(EXC)s.what()')

  # import Eigen3, sva and rbd types
  import_eigen3_types(pg)
  import_sva_types(pg)
  import_rbd_types(pg)

  # build list type
  pg.add_container('std::vector<double>', 'double', 'vector')
  pg.add_container('std::vector<std::vector<double> >', 'std::vector<double>', 'vector')

  # pg
  build_pg(pg)

  with open(sys.argv[1], 'w') as f:
    pg.generate(f)

