#ifndef CLASS_MEMBER_ACCESSOR
#define CLASS_MEMBER_ACCESSOR

#include <cstddef>
#include <iostream>
#include <stdexcept>
#include <string>
#include <typeindex>
#include <unordered_map>

template <class clType> class MemberAccessor {
public:
  // Register a member with its name and type
  template <typename T> static void registerMember(const std::string &name, T clType::*memberPtr) {
    size_t offset         = reinterpret_cast<size_t>(&(((clType *)nullptr)->*memberPtr));
    getMemberInfo()[name] = MemberInfo{typeid(T), offset};
  }

  // Get the offset of a member by name
  static size_t getOffset(const std::string &name) {
    if (getMemberInfo().count(name) == 0) {
      std::cout << "Member not registered: " << name << std::endl;
      return 0;
    }
    return getMemberInfo()[name].offset;
  }

  // Get the type of a member by name
  static std::type_index getType(const std::string &name) {
    if (getMemberInfo().count(name) == 0) {
      std::cout << "Member not registered: " << name << std::endl;
      return typeid(MemberAccessor);
    }
    return getMemberInfo()[name].type;
  }

  // Static helper: Change the value at a given offset for a given instance
  template <typename T> static void setAtOffset(size_t offset, clType &instance, T value) {
    T *memberPtr = reinterpret_cast<T *>(reinterpret_cast<uint8_t *>(&instance) + offset);
    *memberPtr   = value;
  }

  // Static helper: Change the value at a given offset for a given instance
  template <typename T> static void set(const std::string &name, clType &instance, T value) {
    size_t offset = getOffset(name);
    T *memberPtr  = reinterpret_cast<T *>(reinterpret_cast<uint8_t *>(&instance) + offset);
    *memberPtr    = value;
  }

  // Static helper: Get the value at a given offset for a given instance
  template <typename T> static T getAtOffset(size_t offset, const clType &instance) {
    const T *memberPtr = reinterpret_cast<const T *>(reinterpret_cast<const uint8_t *>(&instance) + offset);
    return *memberPtr;
  }

  template <typename T> static T get(const std::string &name, const clType &instance) {
    size_t offset      = getOffset(name);
    const T *memberPtr = reinterpret_cast<const T *>(reinterpret_cast<const uint8_t *>(&instance) + offset);
    return *memberPtr;
  }

private:
  struct MemberInfo {
    std::type_index type;
    size_t offset;
    // Default constructor
    MemberInfo() : type(typeid(void)), offset(0) {}
    // Parameterized constructor
    MemberInfo(std::type_index t, size_t o) : type(t), offset(o) {}
  };

  // Static storage for member info
  static std::unordered_map<std::string, MemberInfo> &getMemberInfo() {
    static std::unordered_map<std::string, MemberInfo> memberInfo;
    return memberInfo;
  }
};

#endif

#if 1

#include <iostream>
#include <vector>

struct Particle {
  double normalStiffness{1e6}, tangentialStiffness{1e6};
  int myInt{42};
  double normalViscDampingRate{0.98}, friction{0.0};
  double rollingFriction{0.123}, adhesion{20.4}, GcGlue{0.0};
  bool myBool{true};
};

int main() {
  // Register members
  MemberAccessor<Particle>::registerMember("friction", &Particle::friction);
  MemberAccessor<Particle>::registerMember("rollingFriction", &Particle::rollingFriction);
  MemberAccessor<Particle>::registerMember("adhesion", &Particle::adhesion);
  MemberAccessor<Particle>::registerMember("myBool", &Particle::myBool);
  MemberAccessor<Particle>::registerMember("myInt", &Particle::myInt);

  Particle p;
  std::vector<Particle> particles(10);

  // Get the offset of 'friction'
  size_t frictionOffset = MemberAccessor<Particle>::getOffset("friction");

  // Set 'friction' for all particles
  for (auto &particle : particles) { MemberAccessor<Particle>::setAtOffset<double>(frictionOffset, particle, 0.5); }

  // Set 'friction' for a single particle
  MemberAccessor<Particle>::setAtOffset<double>(frictionOffset, p, 0.75);

  // Get 'friction' for a single particle
  double friction = MemberAccessor<Particle>::getAtOffset<double>(frictionOffset, p);
  std::cout << "Friction: " << friction << std::endl;

  size_t rfrictionOffset = MemberAccessor<Particle>::getOffset("rollingFriction");
  MemberAccessor<Particle>::setAtOffset<double>(rfrictionOffset, p, 1.65);
  std::cout << "rollingFriction: " << MemberAccessor<Particle>::getAtOffset<double>(rfrictionOffset, p) << std::endl;

  std::cout << "adhesion: " << MemberAccessor<Particle>::get<double>("adhesion", p) << std::endl;
  MemberAccessor<Particle>::set<double>("adhesion", p, -0.32);
  std::cout << "adhesion: " << p.adhesion << std::endl;

  std::cout << "bool: " << MemberAccessor<Particle>::get<bool>("myBool", p) << std::endl;

  std::cout << "int: " << MemberAccessor<Particle>::get<int>("myInt", p) << std::endl;

  std::string who = "myInt";

  std::type_index theType = MemberAccessor<Particle>::getType(who);
  if (theType == typeid(bool)) {
    std::cout << "?? " << MemberAccessor<Particle>::get<bool>(who, p) << std::endl;
  } else if (theType == typeid(int)) {
    std::cout << "?? " << MemberAccessor<Particle>::get<int>(who, p) << std::endl;
  }

  return 0;
}

#endif
