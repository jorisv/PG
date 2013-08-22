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

  fixedPositionContact = pg.add_struct('FixedPositionContact')
  fixedOrientationContact = pg.add_struct('FixedOrientationContact')
  forceContact = pg.add_struct('ForceContact')

  # build list type
  pg.add_container('std::vector<pg::FixedPositionContact>', 'pg::FixedPositionContact', 'vector')
  pg.add_container('std::vector<pg::FixedOrientationContact>', 'pg::FixedOrientationContact', 'vector')
  pg.add_container('std::vector<pg::ForceContact>', 'pg::ForceContact', 'vector')

  # PostureGenerator
  pgSolver.add_constructor([param('const rbd::MultiBody&', 'mb')])

  pgSolver.add_method('fixedPositionContacts', None, [param('std::vector<pg::FixedPositionContact>', 'contacts')])
  pgSolver.add_method('fixedOrientationContacts', None, [param('std::vector<pg::FixedOrientationContact>', 'contacts')])
  pgSolver.add_method('forceContacts', None, [param('std::vector<pg::ForceContact>', 'contacts')])
  pgSolver.add_method('qBounds', None, [param('std::vector<std::vector<double> >', 'lq'),
                                        param('std::vector<std::vector<double> >', 'uq')])


  # Don't change the order. We must try to convert in int before convert in double
  # because double -> int fail but int -> double succeed (so int are read as double).
  pgSolver.add_method('param', None, [param('const std::string&', 'name'), param('const std::string&', 'value')])
  pgSolver.add_method('param', None, [param('const std::string&', 'name'), param('int', 'value')])
  pgSolver.add_method('param', None, [param('const std::string&', 'name'), param('double', 'value')])

  pgSolver.add_method('run', retval('bool'), [param('std::vector<std::vector<double> >', 'q')])
  pgSolver.add_method('q', retval('std::vector<std::vector<double> >'), [])

  # FixedPositionContact
  fixedPositionContact.add_constructor([])

  fixedPositionContact.add_instance_attribute('bodyId', 'int')
  fixedPositionContact.add_instance_attribute('target', 'Eigen::Vector3d')
  fixedPositionContact.add_instance_attribute('surfaceFrame', 'sva::PTransformd')

  # FixedOrientationContact
  fixedOrientationContact.add_constructor([])

  fixedOrientationContact.add_instance_attribute('bodyId', 'int')
  fixedOrientationContact.add_instance_attribute('target', 'Eigen::Matrix3d')
  fixedOrientationContact.add_instance_attribute('surfaceFrame', 'sva::PTransformd')

  # ForceContact
  forceContact.add_constructor([])

  forceContact.add_instance_attribute('bodyId', 'int')
  forceContact.add_instance_attribute('points', 'std::vector<sva::PTransformd>')



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
  pg.add_container('std::vector<sva::PTransformd>', 'sva::PTransformd', 'vector')

  # pg
  build_pg(pg)

  with open(sys.argv[1], 'w') as f:
    pg.generate(f)

