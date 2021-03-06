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


def import_sch_types(mod):
  mod.add_class('S_Object', foreign_cpp_namespace='sch', import_from_module='sch')
  mod.add_class('CD_Pair', foreign_cpp_namespace='sch', import_from_module='sch')


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
  mod.add_class('Vector2d', foreign_cpp_namespace='Eigen', import_from_module='eigen3')
  mod.add_class('Vector6d', foreign_cpp_namespace='Eigen', import_from_module='eigen3')

  mod.add_class('Matrix3d', foreign_cpp_namespace='Eigen', import_from_module='eigen3')
  mod.add_class('Matrix6d', foreign_cpp_namespace='Eigen', import_from_module='eigen3')

  mod.add_class('MatrixXd', foreign_cpp_namespace='Eigen', import_from_module='eigen3')
  mod.add_class('VectorXd', foreign_cpp_namespace='Eigen', import_from_module='eigen3')

  mod.add_class('Quaterniond', foreign_cpp_namespace='Eigen', import_from_module='eigen3')



def build_pg(pg):
  pgSolver = pg.add_class('PostureGenerator', custom_name='PostureGenerator')

  robotConfig = pg.add_struct('RobotConfig')
  fixedPositionContact = pg.add_struct('FixedPositionContact')
  fixedOrientationContact = pg.add_struct('FixedOrientationContact')
  planarContact = pg.add_struct('PlanarContact')
  ellipseContact = pg.add_struct('EllipseContact')
  gripperContact = pg.add_struct('GripperContact')
  forceContact = pg.add_struct('ForceContact')
  envCollision = pg.add_struct('EnvCollision')
  selfCollision = pg.add_struct('SelfCollision')
  bodyPosTarget = pg.add_struct('BodyPositionTarget')
  bodyOriTarget = pg.add_struct('BodyOrientationTarget')
  forceContactMin = pg.add_struct('ForceContactMinimization')
  torqueContactMin = pg.add_struct('TorqueContactMinimization')
  normalForceTarget = pg.add_struct('NormalForceTarget')
  tanForceMin = pg.add_struct('TangentialForceMinimization')
  iterateQuantities = pg.add_struct('IterateQuantities')
  ellipseResult = pg.add_struct('EllipseResult')
  comHalfSpace = pg.add_struct('CoMHalfSpace')
  runConfig = pg.add_struct('RunConfig')
  bodyLink = pg.add_struct('BodyLink')
  robotLink = pg.add_struct('RobotLink')
  cylindricalContact = pg.add_struct('CylindricalContact')

  # build list type
  pg.add_container('std::vector<pg::FixedPositionContact>', 'pg::FixedPositionContact', 'vector')
  pg.add_container('std::vector<pg::FixedOrientationContact>', 'pg::FixedOrientationContact', 'vector')
  pg.add_container('std::vector<pg::PlanarContact>', 'pg::PlanarContact', 'vector')
  pg.add_container('std::vector<pg::EllipseContact>', 'pg::EllipseContact', 'vector')
  pg.add_container('std::vector<pg::GripperContact>', 'pg::GripperContact', 'vector')
  pg.add_container('std::vector<pg::ForceContact>', 'pg::ForceContact', 'vector')
  pg.add_container('std::vector<pg::EnvCollision>', 'pg::EnvCollision', 'vector')
  pg.add_container('std::vector<pg::SelfCollision>', 'pg::SelfCollision', 'vector')
  pg.add_container('std::vector<pg::CoMHalfSpace>', 'pg::CoMHalfSpace', 'vector')
  pg.add_container('std::vector<pg::BodyPositionTarget>', 'pg::BodyPositionTarget', 'vector')
  pg.add_container('std::vector<pg::BodyOrientationTarget>', 'pg::BodyOrientationTarget', 'vector')
  pg.add_container('std::vector<pg::ForceContactMinimization>', 'pg::ForceContactMinimization', 'vector')
  pg.add_container('std::vector<pg::TorqueContactMinimization>', 'pg::TorqueContactMinimization', 'vector')
  pg.add_container('std::vector<pg::NormalForceTarget>', 'pg::NormalForceTarget', 'vector')
  pg.add_container('std::vector<pg::TangentialForceMinimization>', 'pg::TangentialForceMinimization', 'vector')
  pg.add_container('std::vector<Eigen::Vector3d>', 'Eigen::Vector3d', 'vector')
  pg.add_container('std::vector<Eigen::VectorXd>', 'Eigen::VectorXd', 'vector')
  pg.add_container('std::vector<std::vector<Eigen::VectorXd> >', 'std::vector<Eigen::VectorXd>', 'vector')
  pg.add_container('std::vector<pg::EllipseResult>', 'pg::EllipseResult', 'vector')
  pg.add_container('std::vector<pg::RobotConfig>', 'pg::RobotConfig', 'vector')
  pg.add_container('std::vector<pg::RunConfig>', 'pg::RunConfig', 'vector')
  pg.add_container('std::vector<pg::BodyLink>', 'pg::BodyLink', 'vector')
  pg.add_container('std::vector<pg::RobotLink>', 'pg::RobotLink', 'vector')
  pg.add_container('std::vector<pg::CylindricalContact>', 'pg::CylindricalContact', 'vector')

  # PostureGenerator
  pgSolver.add_constructor([])

  pgSolver.add_method('robotConfigs', None, [param('std::vector<pg::RobotConfig>', 'rc'), param('const Eigen::Vector3d&', 'gravity')])
  pgSolver.add_method('robotConfigs', retval('std::vector<pg::RobotConfig>'), [], is_const=True)

  pgSolver.add_method('robotLinks', None, [param('std::vector<pg::RobotLink>', 'rl')])
  pgSolver.add_method('robotLinks', retval('std::vector<pg::RobotLink>'), [], is_const=True)

  # Don't change the order. We must try to convert in int before convert in double
  # because double -> int fail but int -> double succeed (so int are read as double).
  pgSolver.add_method('param', None, [param('const std::string&', 'name'), param('const std::string&', 'value')])
  pgSolver.add_method('param', None, [param('const std::string&', 'name'), param('int', 'value')])
  pgSolver.add_method('param', None, [param('const std::string&', 'name'), param('double', 'value')])

  pgSolver.add_method('run', retval('bool'), [param('const std::vector<pg::RunConfig>&', 'config')])

  pgSolver.add_method('q', retval('std::vector<std::vector<double> >'), [], is_const=True)
  pgSolver.add_method('forces', retval('std::vector<sva::ForceVecd>'), [], is_const=True)
  pgSolver.add_method('torque', retval('std::vector<std::vector<double> >'), [], is_const=True)
  pgSolver.add_method('ellipses', retval('std::vector<pg::EllipseResult>'), [], is_const=True)

  pgSolver.add_method('q', retval('std::vector<std::vector<double> >'),
                      [param('int', 'robot')], is_const=True)
  pgSolver.add_method('forces', retval('std::vector<sva::ForceVecd>'),
                      [param('int', 'robot')], is_const=True)
  pgSolver.add_method('torque', retval('std::vector<std::vector<double> >'),
                      [param('int', 'robot')])
  pgSolver.add_method('ellipses', retval('std::vector<pg::EllipseResult>'),
                      [param('int', 'robot')], is_const=True)

  pgSolver.add_method('nrIters', retval('int'), [], is_const=True)
  pgSolver.add_method('qIter', retval('std::vector<std::vector<double> >'),
                      [param('int', 'iter')], throw=[out_ex], is_const=True)
  pgSolver.add_method('forcesIter', retval('std::vector<sva::ForceVecd>'),
                      [param('int', 'iter')], throw=[out_ex], is_const=True)
  pgSolver.add_method('torqueIter', retval('std::vector<std::vector<double> >'),
                      [param('int', 'iter')], throw=[out_ex])
  pgSolver.add_method('ellipsesIter', retval('std::vector<pg::EllipseResult>'),
                      [param('int', 'iter')], throw=[out_ex], is_const=True)
  pgSolver.add_method('qIter', retval('std::vector<std::vector<double> >'),
                      [param('int', 'robot'), param('int', 'iter')], throw=[out_ex], is_const=True)
  pgSolver.add_method('forcesIter', retval('std::vector<sva::ForceVecd>'),
                      [param('int', 'robot'), param('int', 'iter')], throw=[out_ex], is_const=True)
  pgSolver.add_method('torqueIter', retval('std::vector<std::vector<double> >'),
                      [param('int', 'robot'), param('int', 'iter')], throw=[out_ex])
  pgSolver.add_method('ellipsesIter', retval('std::vector<pg::EllipseResult>'),
                      [param('int', 'robot'), param('int', 'iter')], throw=[out_ex], is_const=True)
  pgSolver.add_method('quantitiesIter', retval('pg::IterateQuantities'),
                      [param('int', 'iter')], throw=[out_ex], is_const=True)

  # FixedPositionContact
  fixedPositionContact.add_constructor([])
  fixedPositionContact.add_constructor([param('int', 'bodyId'),
                                        param('const Eigen::Vector3d&', 't'),
                                        param('const sva::PTransformd&', 'sf')])

  fixedPositionContact.add_instance_attribute('bodyId', 'int')
  fixedPositionContact.add_instance_attribute('target', 'Eigen::Vector3d')
  fixedPositionContact.add_instance_attribute('surfaceFrame', 'sva::PTransformd')

  # FixedOrientationContact
  fixedOrientationContact.add_constructor([])
  fixedOrientationContact.add_constructor([param('int', 'bodyId'),
                                           param('const Eigen::Matrix3d&', 't'),
                                           param('const sva::PTransformd&', 'sf')])

  fixedOrientationContact.add_instance_attribute('bodyId', 'int')
  fixedOrientationContact.add_instance_attribute('target', 'Eigen::Matrix3d')
  fixedOrientationContact.add_instance_attribute('surfaceFrame', 'sva::PTransformd')

  # PlanarContact
  planarContact.add_constructor([])
  planarContact.add_constructor([param('int', 'bodyId'),
                                 param('const sva::PTransformd&', 'tf'), param('std::vector<Eigen::Vector2d>', 'tp'),
                                 param('const sva::PTransformd&', 'sf'), param('std::vector<Eigen::Vector2d>', 'sp')])

  planarContact.add_instance_attribute('bodyId', 'int')
  planarContact.add_instance_attribute('targetFrame', 'sva::PTransformd')
  planarContact.add_instance_attribute('targetPoints', 'std::vector<Eigen::Vector2d>')
  planarContact.add_instance_attribute('surfaceFrame', 'sva::PTransformd')
  planarContact.add_instance_attribute('surfacePoints', 'std::vector<Eigen::Vector2d>')

  # EllipseContact
  ellipseContact.add_constructor([])
  ellipseContact.add_constructor([param('int', 'bodyId'),
                                  param('double', 'radiusMin1'),
                                  param('double', 'radiusMin2'),
                                  param('const sva::PTransformd&', 'tf'),
                                  param('std::vector<Eigen::Vector2d>', 'tp'),
                                  param('const sva::PTransformd&', 'sf'),
                                  param('std::vector<Eigen::Vector2d>', 'sp')])
  ellipseContact.add_constructor([param('int', 'bodyId'),
                                  param('double', 'radiusMin1'),
                                  param('const sva::PTransformd&', 'tf'),
                                  param('std::vector<Eigen::Vector2d>', 'tp'),
                                  param('const sva::PTransformd&', 'sf'),
                                  param('std::vector<Eigen::Vector2d>', 'sp')])

  ellipseContact.add_instance_attribute('bodyId', 'int')
  ellipseContact.add_instance_attribute('radiusMin1', 'double')
  ellipseContact.add_instance_attribute('radiusMin2', 'double')
  ellipseContact.add_instance_attribute('targetFrame', 'sva::PTransformd')
  ellipseContact.add_instance_attribute('targetPoints', 'std::vector<Eigen::Vector2d>')
  ellipseContact.add_instance_attribute('surfaceFrame', 'sva::PTransformd')
  ellipseContact.add_instance_attribute('surfacePoints', 'std::vector<Eigen::Vector2d>')

  # GripperContact
  gripperContact.add_constructor([])
  gripperContact.add_constructor([param('int', 'bodyId'),
                                 param('const sva::PTransformd&', 'tf'), param('std::vector<Eigen::Vector2d>', 'tp'),
                                 param('const sva::PTransformd&', 'sf'), param('std::vector<Eigen::Vector2d>', 'sp')])

  gripperContact.add_instance_attribute('bodyId', 'int')
  gripperContact.add_instance_attribute('targetFrame', 'sva::PTransformd')
  gripperContact.add_instance_attribute('targetPoints', 'std::vector<Eigen::Vector2d>')
  gripperContact.add_instance_attribute('surfaceFrame', 'sva::PTransformd')
  gripperContact.add_instance_attribute('surfacePoints', 'std::vector<Eigen::Vector2d>')

  # ForceContact
  forceContact.add_constructor([])
  forceContact.add_constructor([param('int', 'bodyId'),
                                param('std::vector<sva::PTransformd>', 'points'),
                                param('double', 'mu')])

  forceContact.add_instance_attribute('bodyId', 'int')
  forceContact.add_instance_attribute('points', 'std::vector<sva::PTransformd>')
  forceContact.add_instance_attribute('mu', 'double')

  # EnvCollision
  envCollision.add_constructor([])
  envCollision.add_constructor([param('int', 'bodyId'),
                                param('sch::S_Object*', 'bodyHull', transfer_ownership=False),
                                param('const sva::PTransformd&', 'bodyT'),
                                param('sch::S_Object*', 'envHull', transfer_ownership=False),
                                param('double', 'minDist')])

  envCollision.add_instance_attribute('bodyId', 'int')
  # pybindgen have some issue with ptr return
  # envCollision.add_instance_attribute('bodyHull', retval('sch::S_Object*',caller_owns_return=False))
  envCollision.add_instance_attribute('bodyT', 'sva::PTransformd')
  # envCollision.add_instance_attribute('envHull', retval('sch::S_Object*',caller_owns_return=False))
  envCollision.add_instance_attribute('minDist', 'double')

  # SelfCollision
  selfCollision.add_constructor([])
  selfCollision.add_constructor([param('int', 'body1Id'),
                                 param('sch::S_Object*', 'body1Hull', transfer_ownership=False),
                                 param('const sva::PTransformd&', 'body1T'),
                                 param('int', 'body2Id'),
                                 param('sch::S_Object*', 'body2Hull', transfer_ownership=False),
                                 param('const sva::PTransformd&', 'body2T'),
                                 param('double', 'minDist')])

  selfCollision.add_instance_attribute('body1Id', 'int')
  # pybindgen have some issue with ptr return
  # selfCollision.add_instance_attribute('body2Hull', retval('sch::S_Object*',caller_owns_return=False))
  selfCollision.add_instance_attribute('body1T', 'sva::PTransformd')
  selfCollision.add_instance_attribute('body2Id', 'int')
  # selfCollision.add_instance_attribute('body1Hull', retval('sch::S_Object*',caller_owns_return=False))
  selfCollision.add_instance_attribute('body2T', 'sva::PTransformd')
  selfCollision.add_instance_attribute('minDist', 'double')

  # CoMHalfSpace
  comHalfSpace.add_constructor([])
  comHalfSpace.add_constructor([param('const std::vector<Eigen::Vector3d>&', 'O'),
                                param('const std::vector<Eigen::Vector3d>&', 'n')])
  comHalfSpace.add_instance_attribute('origins', 'std::vector<Eigen::Vector3d>')
  comHalfSpace.add_instance_attribute('normals', 'std::vector<Eigen::Vector3d>')

  # BodyPositionTarget
  bodyPosTarget.add_constructor([])
  bodyPosTarget.add_constructor([param('int', 'bodyId'),
                                 param('const Eigen::Vector3d&', 'target'),
                                 param('double', 'scale')])

  bodyPosTarget.add_instance_attribute('bodyId', 'int')
  bodyPosTarget.add_instance_attribute('target', 'Eigen::Vector3d')
  bodyPosTarget.add_instance_attribute('scale', 'double')

  # BodyOrientationTarget
  bodyOriTarget.add_constructor([])
  bodyOriTarget.add_constructor([param('int', 'bodyId'),
                                 param('const Eigen::Matrix3d&', 'target'),
                                 param('double', 'scale')])

  bodyOriTarget.add_instance_attribute('bodyId', 'int')
  bodyOriTarget.add_instance_attribute('target', 'Eigen::Matrix3d')
  bodyOriTarget.add_instance_attribute('scale', 'double')

  # ForceContactMinimization
  forceContactMin.add_constructor([])
  forceContactMin.add_constructor([param('int', 'bodyId'),
                                   param('double', 'scale')])

  forceContactMin.add_instance_attribute('bodyId', 'int')
  forceContactMin.add_instance_attribute('scale', 'double')

  # TorqueContactMinimization
  torqueContactMin.add_constructor([])
  torqueContactMin.add_constructor([param('int', 'bodyId'),
                                    param('const Eigen::Vector3d&', 'origin'),
                                    param('const Eigen::Vector3d&', 'axis'),
                                    param('double', 'scale')])

  torqueContactMin.add_instance_attribute('bodyId', 'int')
  torqueContactMin.add_instance_attribute('origin', 'Eigen::Vector3d')
  torqueContactMin.add_instance_attribute('axis', 'Eigen::Vector3d')
  torqueContactMin.add_instance_attribute('scale', 'double')

  # NormalForceTarget
  normalForceTarget.add_constructor([])
  normalForceTarget.add_constructor([param('int', 'bodyId'),
                                     param('double', 'target'),
                                     param('double', 'scale')])

  normalForceTarget.add_instance_attribute('bodyId', 'int')
  normalForceTarget.add_instance_attribute('target', 'double')
  normalForceTarget.add_instance_attribute('scale', 'double')

  # TangentialForceMinimization
  tanForceMin.add_constructor([])
  tanForceMin.add_constructor([param('int', 'bodyId'),
                                     param('double', 'scale')])

  tanForceMin.add_instance_attribute('bodyId', 'int')
  tanForceMin.add_instance_attribute('scale', 'double')

  # RobotConfig
  robotConfig.add_constructor([param('const rbd::MultiBody&', 'mb')])
  robotConfig.add_instance_attribute('fixedPosContacts', 'std::vector<pg::FixedPositionContact>')
  robotConfig.add_instance_attribute('fixedOriContacts', 'std::vector<pg::FixedOrientationContact>')
  robotConfig.add_instance_attribute('planarContacts', 'std::vector<pg::PlanarContact>')
  robotConfig.add_instance_attribute('ellipseContacts', 'std::vector<pg::EllipseContact>')
  robotConfig.add_instance_attribute('gripperContacts', 'std::vector<pg::GripperContact>')
  robotConfig.add_instance_attribute('cylindricalContacts', 'std::vector<pg::CylindricalContact>')
  robotConfig.add_instance_attribute('forceContacts', 'std::vector<pg::ForceContact>')
  robotConfig.add_instance_attribute('envCollisions', 'std::vector<pg::EnvCollision>')
  robotConfig.add_instance_attribute('selfCollisions', 'std::vector<pg::SelfCollision>')
  robotConfig.add_instance_attribute('comHalfSpaces', 'std::vector<pg::CoMHalfSpace>')
  robotConfig.add_instance_attribute('ql', 'std::vector<std::vector<double> >')
  robotConfig.add_instance_attribute('qu', 'std::vector<std::vector<double> >')
  robotConfig.add_instance_attribute('tl', 'std::vector<std::vector<double> >')
  robotConfig.add_instance_attribute('tu', 'std::vector<std::vector<double> >')
  robotConfig.add_instance_attribute('tlPoly', 'std::vector<std::vector<Eigen::VectorXd> >')
  robotConfig.add_instance_attribute('tuPoly', 'std::vector<std::vector<Eigen::VectorXd> >')

  robotConfig.add_instance_attribute('postureScale', 'double')
  robotConfig.add_instance_attribute('torqueScale', 'double')
  robotConfig.add_instance_attribute('forceScale', 'double')
  robotConfig.add_instance_attribute('ellipseCostScale', 'double')
  robotConfig.add_instance_attribute('bodyPosTargets', 'std::vector<pg::BodyPositionTarget>')
  robotConfig.add_instance_attribute('bodyOriTargets', 'std::vector<pg::BodyOrientationTarget>')
  robotConfig.add_instance_attribute('forceContactsMin', 'std::vector<pg::ForceContactMinimization>')
  robotConfig.add_instance_attribute('torqueContactsMin', 'std::vector<pg::TorqueContactMinimization>')
  robotConfig.add_instance_attribute('normalForceTargets', 'std::vector<pg::NormalForceTarget>')
  robotConfig.add_instance_attribute('tanForceMin', 'std::vector<pg::TangentialForceMinimization>')

  # RunConfig
  runConfig.add_constructor([])
  runConfig.add_constructor([param('std::vector<std::vector<double> >', 'initQ'),
                             param('std::vector<sva::ForceVecd>', 'initForces'),
                             param('std::vector<std::vector<double> >', 'targetQ')])
  runConfig.add_instance_attribute('initQ', 'std::vector<std::vector<double> >')
  runConfig.add_instance_attribute('initForces', 'std::vector<sva::ForceVecd>')
  runConfig.add_instance_attribute('targetQ', 'std::vector<std::vector<double> >')

  # BodyLink
  bodyLink.add_constructor([])
  bodyLink.add_constructor([param('int', 'bodyId'),
                            param('const sva::PTransformd&', 'body1T'),
                            param('const sva::PTransformd&', 'body2T')])
  bodyLink.add_instance_attribute('bodyId', 'int')
  bodyLink.add_instance_attribute('body1T', 'sva::PTransformd')
  bodyLink.add_instance_attribute('body2T', 'sva::PTransformd')


  # RobotLink
  robotLink.add_constructor([])
  robotLink.add_constructor([param('int', 'robot1Index'),
                             param('int', 'robot2Index'),
                             param('std::vector<pg::BodyLink>', 'linkedBodies')])
  robotLink.add_instance_attribute('robot1Index', 'int')
  robotLink.add_instance_attribute('robot2Index', 'int')
  robotLink.add_instance_attribute('linkedBodies', 'std::vector<pg::BodyLink>')

  # CylindricalContact
  cylindricalContact.add_constructor([])
  cylindricalContact.add_constructor([param('int', 'bodyId'),
                                      param('double', 'targetRadius'),
                                      param('double', 'targetWidth'),
                                      param('sva::PTransformd', 'targetFrame'),
                                      param('sva::PTransformd', 'surfaceFrame')])
  cylindricalContact.add_instance_attribute('bodyId', 'int')
  cylindricalContact.add_instance_attribute('targetRadius', 'double')
  cylindricalContact.add_instance_attribute('targetWidth', 'double')
  cylindricalContact.add_instance_attribute('targetFrame', 'sva::PTransformd')
  cylindricalContact.add_instance_attribute('surfaceFrame', 'sva::PTransformd')

  # IterateQuantities
  iterateQuantities.add_instance_attribute('obj', 'double')
  iterateQuantities.add_instance_attribute('constr_viol', 'double')

  # EllipseResult
  ellipseResult.add_instance_attribute('bodyIndex', 'int')
  ellipseResult.add_instance_attribute('x', 'double')
  ellipseResult.add_instance_attribute('y', 'double')
  ellipseResult.add_instance_attribute('theta', 'double')
  ellipseResult.add_instance_attribute('r1', 'double')
  ellipseResult.add_instance_attribute('r2', 'double')


if __name__ == '__main__':
  if len(sys.argv) < 2:
    sys.exit(1)

  pg = Module('_pg', cpp_namespace='::pg')
  pg.add_include('<PostureGenerator.h>')

  pg.add_include('<sch/S_Object/S_Object.h>')
  pg.add_include('<sch/CD/CD_Pair.h>')

  dom_ex = pg.add_exception('std::domain_error', foreign_cpp_namespace=' ',
                               message_rvalue='%(EXC)s.what()')
  out_ex = pg.add_exception('std::out_of_range', foreign_cpp_namespace=' ',
                               message_rvalue='%(EXC)s.what()')

  # import Eigen3, sva and rbd types
  import_eigen3_types(pg)
  import_sva_types(pg)
  import_rbd_types(pg)
  import_sch_types(pg)

  # build list type
  pg.add_container('std::vector<double>', 'double', 'vector')
  pg.add_container('std::vector<std::vector<double> >', 'std::vector<double>', 'vector')
  pg.add_container('std::vector<sva::PTransformd>', 'sva::PTransformd', 'vector')
  pg.add_container('std::vector<sva::ForceVecd>', 'sva::ForceVecd', 'vector')
  pg.add_container('std::vector<Eigen::Vector2d>', 'Eigen::Vector2d', 'vector')

  # pg
  build_pg(pg)

  with open(sys.argv[1], 'w') as f:
    pg.generate(f)

